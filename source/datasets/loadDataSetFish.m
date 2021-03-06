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
        findCurve( ind(kk), [loadDir, sourceFile], ...
                   [splineDir, splineFile], splineData, ...
                   constSpeed, doPlot );
    end
    
    [d, ~] = loadCurve(splineFile, 'workdir', splineDir);
    dList{kk} = d;
end
          
end

function findCurve( fishInd, sourceFile, splineFile, splineData, ...
                    constSpeed, doPlot )
                
flipLR = [  1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12, ...
           13,  15,  16,  17,  18,  19,  22,  23,  24,  25,  27,  28, ...
           29,  30,  31,  32,  33,  46, 51,  53,  55,  57,  59,  60,  61, ...
           62,  65,  69,  70,  71,  73,  75,  81,  83,  85,  87,  89, ...
           91,  93,  96,  98, 100, 101, 104, 106, 108, 110, 111, 113, ...
          115, 119, 121, 123, 124, 125, 126, 128, 130, 133, 134, 135, ...
          140, 141, 142, 143, 145, 146, 147, 150, 151, 152, 153, 154, ...
          155, 157, 158, 159, 160, 162, 166, 167, 171, 172, 180, 185, ...
          186, 187, 188, 189, 190, 192, 193, 194, 195, 196, 197, 198, ...
          199 ];
                
%% Load curve
C = dlmread(sourceFile, ' ', 1, 0); % delimiter=' '
                                    % row offset=1, column offset=0
                                   
xLen = max(C(:,1)) - min(C(:,1));   % Make fish lie horizontally
yLen = max(C(:,2)) - min(C(:,2));
if yLen > xLen
    C = [C(:,2), -C(:,1)];
end

if ismember(fishInd, flipLR) % Flip horizontally, if necessary
   C(:,1) = -C(:,1); 
end

[~, I] = min(C(:,1));
C = circshift(C, -(I(1)-1), 1); % Start curve at left tip of the fish

d0 = constructSplineApproximation(C, splineData);
[d0, center] = curveCenter(d0, splineData);

if curveArea(d0, splineData) < 0    % Make fish positively oriented
    d0 = flipud(d0);
end

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
plotCurve(d0, splineData, 'lineStyle', 'r');

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
