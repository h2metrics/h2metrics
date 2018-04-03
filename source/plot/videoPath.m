%% videoPath
%
% Creates a video of a path of curves.
%
% Input
%   filename
%       Filename of the video file to be created.
%   dPath
%       The path to be visualized.
%   splineData
%       General information about the splines used.
%
% Optional inputs
%   vidLength = 5 (default)
%       Length of video in seconds.
%   frameRate = 30 (default)
%       Frames per second
%   vidWidth, vidHeight = 560 x 420
%       Size of video in pixels
%   linesStyle, lineWidth
%       Parameters are passed to plotCurve
%
function videoPath(filename, dPath, splineData, varargin)

% Handle optional inputs
p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'vidLength', 5);
addParameter(p, 'frameRate', 30);
addParameter(p, 'vidWidth', 560);
addParameter(p, 'vidHeight', 420);
addParameter(p, 'boundingBox', []);
parse(p, varargin{:});

% Assign optional inputs
vidLength = p.Results.vidLength;
frameRate = p.Results.frameRate;
vidWidth = p.Results.vidWidth;
vidHeight = p.Results.vidHeight;
aspectRatio = vidWidth / vidHeight;
boundingBox = p.Results.boundingBox;

N = vidLength * frameRate; % Number of frames
ptsT = linspace(0, 1, N)';

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

% Create video object
vidObj = VideoWriter(filename);
vidObj.FrameRate = frameRate;

open(vidObj);

close all;
fig = figure('Renderer', 'painters', ...
             'Position', [1, 1, vidWidth, vidHeight]);

plot(xmin-1, ymin-1); % Create empty plot
axis([xmin, xmax, ymin, ymax]);

for jj = 1:N
    d = evalPath(dPath, ptsT(jj), splineData);
    
    % Create plot
    hold on;
    plotCurve(d, splineData, varargin{:});
    hold off;
    
    axis([xmin, xmax, ymin, ymax]);
    
    drawnow;
    pause(0.05);

    % Store the frame
    newFrameOut=getframe; % leaving gcf out crops the frame in the movie.
    writeVideo(vidObj, newFrameOut);
    
    % Clear figure for next frame
    cla(fig);
end

% Output the movie as an avi file
close(vidObj);
close(fig);

end