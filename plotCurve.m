function plotCurve(d0, splineData)

dSpace = splineData.dSpace;
N = splineData.N;
nS = splineData.nS;
knotsS = splineData.knotsS;
quadData = splineData.quadData;

% Plot parameters
noPlotPtsS = 500;
plotPtsS = linspace(0, 2*pi, noPlotPtsS+1);
        
c0 = deBoor(knotsS, nS, d0, plotPtsS, 1, 'periodic', true);
pt0 = deBoor(knotsS, nS, d0, 0, 1, 'periodic', true);

%% Setup plotting
hold on;
axis equal;

%% Do plotting
plot(c0(:, 1), c0(:, 2), 'k-', 'LineWidth', 1);
plot(pt0(1), pt0(2), 'ko', 'LineWidth', 1);
