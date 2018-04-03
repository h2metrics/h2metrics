%% findCurveBoundingBox
%
% Finds the smalles rectangular box to fit all the curves in dList
%
% Input
%   dList
%       List of curves. Can be a cell array or a matrix of dimension
%         [N, 2, noCurves]
%   splineData
%       General information about the splines used.
%
% Output
%   xmin, ymin, xmax, ymax
%       Coordinates of the bounding box
%
function [xmin, ymin, xmax, ymax] = findCurveBoundingBox(dList, splineData)

% Plot parameters
noPlotPts = 100;
plotPts = linspace(0, 2*pi, noPlotPts+1);

% Treat everything as a cell array
if ndims(dList) == 3
    dTmp = {};
    for jj = size(dList,3):-1:1
        dTmp{jj} = dList(:,:,jj);
    end
    dList = dTmp;
elseif ~isa(dList, 'cell')
    dList = {dList}; 
end  

xmin = Inf;
ymin = Inf;
xmax = -Inf;
ymax = -Inf;

for jj = 1:length(dList)
    d = dList{jj};
    
    pts = evalCurve(plotPts, d, splineData);
    xmin = min(xmin, min(pts(:,1)));
    xmax = max(xmax, max(pts(:,1)));
    ymin = min(ymin, min(pts(:,2)));
    ymax = max(ymax, max(pts(:,2)));
end

end