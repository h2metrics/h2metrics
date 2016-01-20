%% loadDataSetHela
% Loads the HeLa cell data of Murphy et al.
%
% For allowed parameters see documentation of loadDataSet.
function dList = loadDataSetHela( splineData, dataDir, varargin )

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

splineDir = [ dataDir, 'splines/hela_murphy/n', ...
              num2str(splineData.nS), '_N', ...
              num2str(splineData.N), constSpeedSuffix, '/'];
loadDir = [ dataDir, 'source/hela_murphy/' ];

listAll = dir([loadDir, 'r*--1---2.dat.png']);
noCurves = length(listAll);

if isempty(ind)
    ind = 1:noCurves;
end

dList = {};
for kk = length(ind):-1:1
    sourceFile = listAll(ind(kk)).name;
    splineFile = ['cell_', num2str(ind(kk), '%02u')];
    
    if reloadData || ~exist([splineDir, splineFile, '.mat'], 'file')
        helaFindCurve( [loadDir, sourceFile], ...
                       [splineDir, splineFile], splineData, ...
                       constSpeed, plotCurve );
    end
    
    [d, ~] = loadCurve(splineFile, 'workdir', splineDir);
    dList{kk} = d;
end
          
end

function helaFindCurve( sourceFile, splineFile, splineData, ...
                        constSpeed, plotCurve )

%% Setup interpolation data
% This spline is used to interpolate the boundary
auxN = 12;
auxnS = 4;
auxSplineData = constructEmptySplineData;
auxSplineData.N = auxN;
auxSplineData.nS = auxnS;
auxSplineData.dSpace = splineData.dSpace;
auxSplineData = constructKnots(auxSplineData);

% This spline data has interpolation parametrs set; don't use for saving
splineDataInterpol = splineData;
splineDataInterpol.noInterpolS = 12 * splineData.N;
splineDataInterpol = constructKnots(splineDataInterpol);
splineDataInterpol = setupQuadData(splineDataInterpol);
quadDataInterpol = splineDataInterpol.quadData;

%% Load curve
% Use Otsu's method of thresholding
I = imread(sourceFile);
eff_luminance = graythresh(I);
% disp(eff_luminance);

BW = im2bw(I, eff_luminance);
[B, ~] = bwboundaries(BW, 4, 'noholes');
lengthComp = cellfun(@(t) size(t, 1), B);
[~, ind] = max(lengthComp);  

% Extract boundary
boundary = B{ind};
boundary = flip(boundary, 2); % Switch x- and y-coordinates

% Do interpolation; interpolate first with auxSplineData, then
% upsample to higher order
boundary = boundary(1:end-1,:); % Last point equals first
d0 = constructSplineApproximation(boundary, auxSplineData);
dPts = deBoor(auxSplineData.knotsS, auxSplineData.nS, d0, ...
    splineDataInterpol.interpolS, 1, 'periodic', true);
d0 = quadDataInterpol.B_interpolS \ dPts;

% Reparametrize to constant speed
if constSpeed
    d0 = curveReparamConstSpeed(d0, splineDataInterpol);
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
plot(boundary(:,1), boundary(:,2) , 'b', 'LineWidth', 2);
hold on;
axis manual;
imshow(I);
plot(boundary(:,1), boundary(:,2) , 'b', 'LineWidth', 2);
hold off;

% Superimpose interpolation over leaf
noPlotPoints = size(boundary, 1);
plotPoints = linspace(0, 2*pi, noPlotPoints);
boundary2 = deBoor( splineData.knotsS, splineData.nS, d0, ...
                    plotPoints, 1, 'periodic', true);
hold on;
plot(boundary2(:,1), boundary2(:,2) , 'r', 'LineWidth', 2);
hold off;

% Save figure
figname = [splineFile, '.jpg'];
% export_fig(figname);
saveTightFigure(handle, figname);
end
