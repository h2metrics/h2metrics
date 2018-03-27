%% plotCurve
%
% Helper function to plot a curve
%
% Input
%   d
%       A curve or a cell array of curves
%   splineData
%       General information about the splines used.
%
% Optional inputs
%   lineStyle = 'k-' (default)
%       lineStyle = 'colour' gives the parametrization plotted along
%           the curve via a colormap
%       Other lineStyle parameters are passed to plot.
%   lineWidth = 1 (default)
%       Passed to plot. With lineStyle='colour', recommended lineWidth=2.
%
function plotCurve(d, splineData, varargin)

% Handle optional inputs
p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'noPts', []);
addParameter(p, 'lineStyle', 'k-');
addParameter(p, 'lineWidth', 1);
parse(p, varargin{:});

% Assign optional inputs
noPts = p.Results.noPts;
if isempty(noPts)
    noPts = 5 * splineData.N;
end
lineStyle = p.Results.lineStyle;
lineWidth = p.Results.lineWidth;

% Plot parameters
plotPtsS = linspace(0, 2*pi, noPts);

% Treat everything as a cell array
if ~isa(d, 'cell')
    d = {d}; 
end        
if ~isa(lineStyle, 'cell')
    lineStyle = {lineStyle};
end

noCurves = length(d);
for ii = 1:noCurves
    if length(lineStyle) < ii
        lS = lineStyle{end};
    else
        lS = lineStyle{ii}; 
    end
    
    c0 = evalCurve(plotPtsS, d{ii}, splineData);
    pt0 = evalCurve(0, d{ii}, splineData);
    
    % Setup plotting
    hold on;
    
    % Do plotting
    if strcmp(lS, 'colour')
        cm=hsv(noPts);
    
        for jj=1:noPts-1
            plot(c0(jj:jj+1,1), c0(jj:jj+1,2), ...
                 'color', cm(jj,:), 'lineWidth', lineWidth);
        end
    else
        plot(c0(:, 1), c0(:, 2), lS);
        plot(pt0(1), pt0(2), 'ko');      
    end
    
    hold off; 
end
end


