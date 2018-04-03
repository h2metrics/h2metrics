%% plotPathRow
%
% Plots snapshots from a path aligned horizontally
%
% Input
%   dPath
%       The path in question
%   splineData
%       General information about the splines used.
%
% Optional inputs
%   noPtsT = 5 (default)
%       How many curves should be plotted
%   blockWidth, blockHeight = 200 (default)
%       Size (in pixels) of each curve
%   marginX, marginY
%       Optional margin (relative to width/height) to be added to the
%       bounding box.
%   lineStyle, lineWidth
%       Parameters passed on to plotCurve
%
function fig = plotPathRow(dPath, splineData, varargin)

% Handle optional inputs
p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'noPtsT', 5);
addParameter(p, 'blockWidth', 200);
addParameter(p, 'blockHeight', 200);
addParameter(p, 'boundingBox', []);
addParameter(p, 'marginX', 0.);
addParameter(p, 'marginY', 0.);
parse(p, varargin{:});

noPtsT = p.Results.noPtsT;
blockWidth = p.Results.blockWidth;
blockHeight = p.Results.blockHeight;
aspectRatio = blockWidth / blockHeight;
boundingBox = p.Results.boundingBox;
marginX = p.Results.marginX;
marginY = p.Results.marginY;

ptsT = linspace(0, 1, noPtsT);

% Find maximum dimensions of all curves
if ~isempty(boundingBox)
    xmin = boundingBox(1);
    ymin = boundingBox(2);
    xmax = boundingBox(3);
    ymax = boundingBox(4);
else
    [xmin,ymin,xmax,ymax] = findCurveBoundingBox(...
        evalPath(dPath, ptsT, splineData), splineData);
end

% Add margin
xmin = xmin - marginX * (xmax - xmin);
xmax = xmax + marginX * (xmax - xmin);
ymin = ymin - marginY * (ymax - ymin);
ymax = ymax + marginY * (ymax - ymin);

% Adjust maximum limits depending on the aspect ratio
if xmax-xmin > aspectRatio*(ymax-ymin)
    delta = 0.5 * ((xmax-xmin)/aspectRatio - (ymax-ymin));
    ymin = ymin - delta;
    ymax = ymax + delta;
else
    delta = 0.5 * (aspectRatio*(ymax-ymin) - (xmax-xmin));
    xmin = xmin - delta;
    xmax = xmax + delta;
end

fig = figure('Renderer', 'painters', ...
             'Position', [1, 1, blockWidth * noPtsT, blockHeight], ...
             'Color', 'white');

plot(xmin-1, ymin-1); % Create empty plot

axis([xmin, xmin + noPtsT * (xmax - xmin), ymin, ymax]);

for jj = 1:noPtsT
    d = evalPath(dPath, ptsT(jj), splineData);
    
    dTmp = d;
    dTmp(:,1) = dTmp(:,1) + (jj-1) * (xmax-xmin);
    
    plotCurve(dTmp, splineData, varargin{:});
end

axis off;
axis equal;

end   