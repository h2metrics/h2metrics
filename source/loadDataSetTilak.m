%% loadDataSet
% Loads the corpus callosum data of Tilak.
%
% Optional arguments
%   'class' = {'ad', 'normal', 'all' (default)}
%       Returns only diseased, healthy or all patients. If class is 'all',
%       then dList will be a cell array with
%           dList{1} = cell array of AD patients
%           dList{2} = cell array of normal patients
%       Also, if indices are specified, then ind has to be a cell array
%       with ind{1} and ind{2} describing the indices of AD and normal
%       patients respectively. Example:
%           loadDataSetTilak( ..., 'class', 'ad', 'ind', [1, 2] )
%           loadDataSetTilak( ..., 'class', 'all', 'ind', {[1, 2], [1]} )
%
function dList = loadDataSetTilak( splineData, dataDir, varargin )

%% Optional arguments
constSpeed = false; % Reparametrize to constant speed
reloadData = false; % Redo extraction
plotCurve = true;
ind = []; % We want only these curves
patientClass = 'all';

ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'constspeed'
                constSpeed = true;
            case 'reloaddata'
                reloadData = true;
            case 'noplot'
                plotCurve = false;
            case 'ind'
                ii = ii + 1;
                ind = varargin{ii};
            case 'class'
                ii = ii + 1;
                patientClass = varargin{ii};
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    end
    ii = ii + 1;
end

%% Prepare directories
if ~isempty(dataDir) && dataDir(end) ~= '/'
    dataDir(end+1) = '/';
end

if constSpeed
    constSpeedSuffix = '_csp';
else
    constSpeedSuffix = '';
end

splineDir = [ dataDir, 'splines/corpus_callosum_tilak/', ...
             'n', num2str(splineData.nS), '_N', ...
              num2str(splineData.N), constSpeedSuffix, '/' ];
loadDir = [ dataDir, 'source/corpus_callosum_tilak/' ];

loadAdDir = [ loadDir, 'OAS1_cc_AD_txt/' ];
splineAdDir = [ splineDir, 'ad/'];
if ~exist(splineAdDir, 'dir')
    mkdir(splineAdDir);
end

loadNormalDir = [ loadDir, 'OAS1_cc_normal_txt/' ];
splineNormalDir = [ splineDir, 'normal/' ];
if ~exist(splineNormalDir, 'dir')
    mkdir(splineNormalDir);
end

switch lower(patientClass)
    case 'ad'
        dList = tilakLoadClass( loadAdDir, splineAdDir, splineData, ...
            ind, constSpeed, plotCurve, reloadData );
    case 'normal'
        dList = tilakLoadClass( loadNormalDir, splineNormalDir, ...
            splineData, ind, constSpeed, plotCurve, reloadData );
    case 'all'
        if isempty(ind)
            ind = {[], []};
        end
        dList{1} = tilakLoadClass( loadAdDir, splineAdDir, splineData, ...
            ind{1}, constSpeed, plotCurve, reloadData );
        dList{2} = tilakLoadClass( loadNormalDir, splineNormalDir, ...
            splineData, ind{2}, constSpeed, plotCurve, reloadData );
end

end

function dList = tilakLoadClass( loadDir, saveDir, splineData, ind, ...
                                 constSpeed, plotCurve, reloadData )

listAll = dir([loadDir, 'OAS1_*.txt']);
noCurves = length(listAll);

if isempty(ind)
    ind = 1:noCurves;
end

dList = {};
for kk = length(ind):-1:1
    sourceFile = listAll(ind(kk)).name;
    splineFile = sourceFile(1:end-4); % No extension for filename
    
    if reloadData || ~exist([saveDir, splineFile, '.mat'], 'file')
        tilakFindCurve( [loadDir, sourceFile], ...
                        [saveDir, splineFile], splineData, ...
                        constSpeed, plotCurve );
    end
    
    [d, ~] = loadCurve(splineFile, 'workdir', saveDir);
    dList{kk} = d;
end

end

function tilakFindCurve( sourceFile, splineFile, splineData, ...
                         constSpeed, plotCurve )

C = dlmread(sourceFile);
C = C(2:end, :); % First row contains only number of data points
d0 = constructSplineApproximation(C, splineData);
[d0, center] = curveCenter(d0, splineData);

% Reparametrize to constant speed
if constSpeed
    d0 = curveReparamConstSpeed(d0, splineData);
end

% Save curve
saveCurve(splineFile, d0, splineData, 'manualName', true);

% Plotting
if ~plotCurve
    return
end

% Superimpose boundary over cell -- curve
handle = figure(2);
handle.Visible = 'off';
plot(C(:,1)-center(1), C(:,2)-center(2) , 'b', 'LineWidth', 2);

% Superimpose interpolation over leaf
noPlotPoints = size(C, 1);
plotPoints = linspace(0, 2*pi, noPlotPoints);
boundary2 = deBoor( splineData.knotsS, splineData.nS, d0, ...
                    plotPoints, 1, 'periodic', true);
pt0 = deBoor( splineData.knotsS, splineData.nS, d0, ...
              0, 1, 'periodic', true );
hold on;
plot(boundary2(:,1), boundary2(:,2) , 'r', 'LineWidth', 2);
plot(pt0(1), pt0(2), 'ko', 'LineWidth', 2);
hold off;

% Save figure
figname = [splineFile, '.jpg'];
export_fig(figname);
% saveTightFigure(handle, figname);

end