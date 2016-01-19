function plotTangentVector(d, v, splineData)

N = splineData.N;
nS = splineData.nS;
knotsS = splineData.knotsS;

% Plot parameters
noPlotPtsS = 500;

plotPtsS = linspace(0, 2*pi, noPlotPtsS+1);

c = deBoor(knotsS, nS, d, plotPtsS, 1, 'periodic', true);
w = deBoor(knotsS, nS, v, plotPtsS, 1, 'periodic', true);

%% Setup plotting
hold on;
axis equal;

%% Start  plotting
plot(c(:, 1), c(:, 2), 'k-', 'LineWidth', 1);
quiver(c(:, 1), c(:, 2), w(:, 1), w(:, 2), 2, 'k-', 'LineWidth', 1);

%% Finish plotting
hold off;