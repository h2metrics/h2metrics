%% plotTangentVector
%
% Helper function to plot a vector field along a curve
%
% Input
%   d
%       The curve
%   v
%       The vector field
%   splineData
%       General information about the splines used.
%
function plotTangentVector(d, v, splineData)

% Plot parameters
noPlotPtsS = 100;

plotPtsS = linspace(0, 2*pi, noPlotPtsS+1);

c = evalCurve(plotPtsS, d, splineData);
w = evalCurve(plotPtsS, v, splineData);

%% Setup plotting
hold on;
axis equal;

%% Start  plotting
plot(c(:, 1), c(:, 2), 'k-', 'LineWidth', 1);
quiver(c(:, 1), c(:, 2), w(:, 1), w(:, 2), 2, 'k-', 'LineWidth', 1);

%% Finish plotting
hold off;