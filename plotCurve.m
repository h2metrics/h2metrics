%% plotCurve
%
% Helper function to plot a curve
%
% Input
%   d0
%       The curve
%   splineData
%       General information about the splines used.
%
% Optional inputs
%   lineStyle = 'k-' (default)
%       lineStyle parameter to be passed to plot.
%
function plotCurve(d0, splineData, lineStyle)

if nargin < 3
    lineStyle = 'k-';
end

nS = splineData.nS;
knotsS = splineData.knotsS;

% Plot parameters
noPlotPtsS = 500;
plotPtsS = linspace(0, 2*pi, noPlotPtsS+1);
        
c0 = deBoor(knotsS, nS, d0, plotPtsS, 1, 'periodic', true);
pt0 = deBoor(knotsS, nS, d0, 0, 1, 'periodic', true);

%% Setup plotting
hold on;
axis equal;

%% Do plotting
plot(c0(:, 1), c0(:, 2), lineStyle, 'LineWidth', 1);
plot(pt0(1), pt0(2), 'ko', 'LineWidth', 1);

hold off;
