%% loadDataSetWings
% Loads the mosquito wings data of
%   F. James Rohlf, James W. Archie, ``A Comparison of Fourier Methods for
%   the Description of Wing Shape in Mosquitoes (Diptera: Culicidae)''. 
%   Syst. Zool., 33(3):302-317, 1984.
%
% For allowed parameters see documentation of loadDataSet.
function dList = loadDataSetWings( splineData, dataDir, varargin )

%% Optional arguments
constSpeed = false;
reloadData = false;
plotCurve = true;
ind = []; % We want only these curves

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

splineDir = [ dataDir, 'splines/mosquito_wings/n', ...
              num2str(splineData.nS), '_N', ...
              num2str(splineData.N), constSpeedSuffix, '/'];
if ~exist(splineDir, 'dir')
    mkdir(splineDir);
end          

%% Load data
dataFile = [dataDir, 'source/mosquito_wings/RohlfArchieWingOutlines.nts'];

dataFileId = fopen(dataFile);

textscan(dataFileId, '%s', 5, 'Delimiter','\n');

noRows = 2537 - 5;
noCols = 10;

dataPts = zeros(noRows * noCols, 1);
for kk = 1:noRows
    if kk == 6
        disp('a');
    end
    C = textscan(dataFileId, '%8c', 10, 'Delimiter', '', 'Whitespace', '');
    pts = C{1};
    for jj = 1:noCols
        dataPts(noCols*(kk-1) + jj) = str2double(pts(jj,:));
    end
end

noCurves = 127 - 1; % Last one seems to be incomplete
dSpace = 2;
noPts = 100;

dataPts = reshape(dataPts(1:dSpace*noPts*noCurves), ...
                  dSpace, noPts, noCurves);
dataPts = permute(dataPts, [2, 1, 3]);

%% Create splines as necessary
if isempty(ind)
    ind = 1:noCurves;
else
    if min(ind(:)) < 1 || max(ind(:)) > noCurves
        error('Invalid index passed to loadDataSetWings.');
    end
end

dList = {};
for kk = length(ind):-1:1
    splineFile = ['wing_', num2str(ind(kk), '%03u')];
    
    if reloadData || ~exist([splineDir, splineFile, '.mat'], 'file')
        wingFindCurve( dataPts(:, :, ind(kk)), ...
                       [splineDir, splineFile], splineData, ...
                       constSpeed, plotCurve );
    end
    
    [d, ~] = loadCurve(splineFile, 'workdir', splineDir);
    dList{kk} = d;
end
          
end

function wingFindCurve( curvePts, splineFile, splineData, ...
                        constSpeed, doPlot )

% Find spline approximation
d0 = constructSplineApproximation(curvePts, splineData);
[d0, center] = curveCenter(d0, splineData);

% Reparametrize to constant speed
if constSpeed
    d0 = curveReparamConstSpeed(d0, splineData);
end

% Save curve
saveCurve(splineFile, d0, splineData, 'manualName', true);

% Plotting
if ~doPlot
    return
end

handle = figure(2);
handle.Visible = 'off';
set(handle, 'defaultLineLineWidth', 2);

% Original curve and spline
plot( curvePts(:,1)-center(1), curvePts(:,2)-center(2) , 'b');  
plotCurve(d0, splineData, 'r');

% Start point of spline
pt0 = deBoor( splineData.knotsS, splineData.nS, d0, ...
              0, 1, 'periodic', true );

hold on;
plot(pt0(1), pt0(2), 'ko');
hold off;
          
% Save figure
figname = [splineFile, '.jpg'];
export_fig(figname);

end
