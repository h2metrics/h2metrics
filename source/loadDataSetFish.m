%% loadDataSetFish
%
% Loads the surrey fish database. Database contains 1100 fish contours each
% contour consisting of 400 to 1600 points. Details can be found in
%
%   Mokhtarian, F. Abbasi S. and Kittler J. ``Robust and Efficient Shape 
%   Indexing through Curvature Scale Space '' in Proceedings of the sixth 
%   British Machine Vision Conference, BMVC'96. Edinburgh, 10-12 September 
%   1996, pp 53-62.
%
%   Mokhtarian, F. Abbasi S. and Kittler J. ``Efficient and Robust 
%   Retrieval by Shape Content through Curvature Scale Space '' in 
%   Proceedings of the First International Workshop on Image Database 
%   and Multimedia Search, Amsterdam, The Netherlands Aug 1996, pp 35-42. 
% 
% For allowed parameters see documentation of loadDataSet.
function dList = loadDataSetFish( splineData, dataDir, varargin )

%% Optional arguments
constSpeed = false;
reloadData = false;
doPlot = true;
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
                doPlot = false;
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

splineDir = [ dataDir, 'splines/surrey_fish/n', ...
              num2str(splineData.nS), '_N', ...
              num2str(splineData.N), constSpeedSuffix, '/'];
if ~exist(splineDir, 'dir')
    mkdir(splineDir);
end          
          
loadDir = [ dataDir, 'source/surrey_fish/' ];

noCurves = 1100; % Fish with index kk has filename
                 % 'kk' + num2str(kk) + '.c'

if isempty(ind)
    ind = 1:noCurves;
end

dList = {};
for kk = length(ind):-1:1
    sourceFile = ['kk', num2str(ind(kk)), '.c'];
    splineFile = ['fish_', num2str(ind(kk), '%04u')];
    
    if reloadData || ~exist([splineDir, splineFile, '.mat'], 'file')
        findCurve( [loadDir, sourceFile], ...
                   [splineDir, splineFile], splineData, ...
                   constSpeed, doPlot );
    end
    
    [d, ~] = loadCurve(splineFile, 'workdir', splineDir);
    dList{kk} = d;
end
          
end

function findCurve( sourceFile, splineFile, splineData, ...
                    constSpeed, doPlot )

%% Load curve
C = dlmread(sourceFile, ' ', 1, 0); % delimiter=' '
                                    % row offset=1, column offset=0
C = fliplr(C); 
C = -C;        % Make fish lie horizontally
d0 = constructSplineApproximation(C, splineData);
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
handle.Visible = 'on';
set(handle, 'defaultLineLineWidth', 2);

% Original curve and spline
plot( C(:,1)-center(1), C(:,2)-center(2) , 'b');  
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
