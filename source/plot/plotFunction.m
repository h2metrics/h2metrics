%% plotFunction
%
% Helper function to plot a function
%
% Input
%   d
%       A curve or a cell array of curves
%   splineData
%       General information about the splines used.
%
% Optional inputs
%   lineSpec = 'k-' (default)
%       lineSpec parameter to be passed to plot.
%
function plotFunction(f, splineData, varargin)

% Handle optional inputs
p = inputParser;
addOptional(p, 'lineSpec', 'k-');
addParameter(p, 'noPts', 300);
parse(p, varargin{:});
r = p.Results;

% Plot parameters
noPlotPtsS = r.noPts;
plotPtsS = linspace(0, 2*pi, noPlotPtsS+1);
lS = r.lineSpec;

% Treat everything as a cell array
if ~isa(f, 'cell')
    f = {f}; 
    lS = {lS};
elseif ~isa(lS, 'cell')
    tmp = lS; lS = cell(length(f), 1);
    [lS{1:length(f)}] = deal(tmp);
end

noCurves = length(f);
for ii = 1:noCurves
    
    f0 = evalCurve(plotPtsS, f{ii}, splineData);
    
    %% Setup plotting
    hold on;
    %axis equal;
    
    %% Do plotting
    plot(plotPtsS, f0, lS{ii});
    
    hold off; 
end
end


