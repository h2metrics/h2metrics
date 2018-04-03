%% loadDataSetXavier
% Loads the heart curves of Xavier
%
% Optional arguments
%   'version' = {'v1' (default), 'v2'}
%       Describes, which version of the data to use.
%   'ind'='noLoops' (only if version='v1')
%       Returns only curves, which do not have selfintersections. Usage
%           loadDataSet( 'xavier_heart', ..., 
%                        'version', 'v1', 'ind', 'noLoops' );
%
function dList = loadDataSetXavier(splineData, dataDir, varargin)

%% Optional arguments
constSpeed = false; % Reparametrize to constant speed
reloadData = false; % Redo extraction
plotCurve = true;
ind = []; % We want only these curves
curveVersion = 'v1';

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
            case 'version'
                ii = ii + 1;
                curveVersion = varargin{ii};
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

splineDir = [ dataDir, 'splines/xavier_heart/', ...
             'n', num2str(splineData.nS), '_N', ...
              num2str(splineData.N), constSpeedSuffix, '/' ];
loadDir = [ dataDir, 'source/xavier_heart/' ];

switch lower(curveVersion)
    case 'v1'
        splineDir = [splineDir, 'v1/'];
        loadDir = [loadDir, 'Curves_v1/'];
    case 'v2'
        splineDir = [splineDir, 'v2/'];
        loadDir = [loadDir, 'Curves_v2/'];
    otherwise
        error('Wrong argument for curve version');
end

if ~exist(splineDir, 'dir')
    mkdir(splineDir);
end

%% Load curves
listAll = dir([loadDir, 'Patient_*']);
noCurves = length(listAll);

if isempty(ind)
    ind = 1:noCurves;
elseif isa(ind, 'char') && strcmpi(ind, 'noloops')
    if ~strcmpi(curveVersion, 'v1')
        error('noLoops works only for v1');
    end
    ind = [     2,      4,  5,  6,      8,  9, 10,     12,  13, ...
           14,     16, 17, 18, 19, 20, 21, 22, 23, 24, 25 ];
end

dList = {};
for kk = length(ind):-1:1
    sourceFile = ['Patient_', num2str(ind(kk))];
    splineFile = sourceFile; % Xavier files have no extension
    
    if reloadData || ~exist([splineDir, splineFile, '.mat'], 'file')
        xavierFindCurve( [loadDir, sourceFile], ...
                         [splineDir, splineFile], splineData, ...
                         constSpeed, plotCurve );
    end
    
    [d, ~] = loadCurve(splineFile, 'workdir', splineDir);
    dList{kk} = d;
end

end

function xavierFindCurve( sourceFile, splineFile, splineData, ...
                          constSpeed, plotCurve )
                     
C3d = dlmread(sourceFile);

%% Transform affine 3d coordinates to 2d coordinates
eVp = [ 1/sqrt(3) sqrt(6)/3 0; 1/sqrt(3) -sqrt(6)/6 1/sqrt(2); ...
        1/sqrt(3) -sqrt(6)/6 -1/sqrt(2)]; %p to e change, orthogonal matrix
pVe = eVp'; %e to b change
%(pVe.v)' = v'*eVp

C = (C3d - ones([size(C3d, 1), 1]) * [1 0 0]) * eVp;
C = C(:,2:3); % The first coordinate is 0
C = C(1:end-1,:); % Last point equals first

%% Construct spline approximation
if splineData.N <= size(C, 1)
    d0 = constructSplineApproximation(C, splineData);
else
    auxSplineData = constructEmptySplineData;
    auxSplineData.N = size(C, 1);
    auxSplineData.nS = splineData.nS;
    auxSplineData.dSpace = 2;
    auxSplineData = constructKnots(auxSplineData);
    quadData = splineData.quadData;
    
    dTemp = constructSplineApproximation(C, auxSplineData);
    dPts = deBoor(auxSplineData.knotsS, auxSplineData.nS, dTemp, ...
        splineData.interpolS, 1, 'periodic', true);
    d0 = quadData.B_interpolS \ dPts;
end

% Reparametrize to constant speed
if constSpeed
    d0 = curveReparamConstSpeed(d0, splineData);
end

%% Save curve
saveCurve(splineFile, d0, splineData, 'manualName', true);

%% Plotting
if ~plotCurve
    return
end

% Superimpose boundary over cell -- curve
handle = figure(2);
handle.Visible = 'off';
plot(C(:,1), C(:,2) , 'bx', 'LineWidth', 2);

% Superimpose interpolation over leaf
noPlotPoints = 100;
plotPoints = linspace(0, 2*pi, noPlotPoints);
boundary2 = deBoor( splineData.knotsS, splineData.nS, d0, ...
                    plotPoints, 1, 'periodic', true);
hold on;
plot(boundary2(:,1), boundary2(:,2) , 'k', 'LineWidth', 2);
hold off;

% Save figure
figname = [splineFile, '.jpg'];
export_fig(figname);
% saveTightFigure(handle, figname);

end