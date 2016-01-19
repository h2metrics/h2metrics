function plotPath(dPath, splineData)

dSpace = splineData.dSpace;
N = splineData.N;
Nt = splineData.Nt;
nS = splineData.nS;
nT = splineData.nT;
knotsS = splineData.knotsS;
knotsT = splineData.knotsT;
quadData = splineData.quadData;

% Plot parameters
noPlotPtsS = 500;
noPlotPtsT = 500;
noSnapshotsT = 5;
noLinesS = 20;

plotPtsS = linspace(0, 2*pi, noPlotPtsS+1);
plotPtsT = linspace(0, 1, noPlotPtsT);
linePtsS = linspace(0, 2*pi, noLinesS+1);
linePtsS = linePtsS(1:end-1);
snapshotPtsT = linspace(0, 1, noSnapshotsT + 2);
snapshotPtsT = snapshotPtsT(2:end-1);
        
d0 = dPath(1:N, :);
d1 = dPath(end-N+1:end, :);

c0 = deBoor(knotsS, nS, d0, plotPtsS, 1, 'periodic', true);
c1 = deBoor(knotsS, nS, d1, plotPtsS, 1, 'periodic', true);

%% Setup plotting
hold on;
axis equal;
% axis tight;
% axis off;

%% Start  plotting
plot(c0(:, 1), c0(:, 2), 'k-', 'LineWidth', 1);
plot(c1(:, 1), c1(:, 2), 'k-', 'LineWidth', 1);

% Particle paths across time
for ii = 1:noLinesS
    for ll = Nt:-1:1
        lineT(ll, :) = deBoor(knotsS, nS, ...
            dPath((ll-1)*N+1:ll*N,:), linePtsS(ii), 1, ...
            'periodic', true);
    end
    b = deBoor(knotsT, nT, lineT, ...
               plotPtsT, 1, 'periodic', false);
    plot(b(:, 1), b(:, 2), 'k-', 'LineWidth', 1);
end

% Snapshots of curve
curveS = [];
for ii = 1:noSnapshotsT
    for ll = N:-1:1
        curveS(ll, :) = deBoor(knotsT, nT, ...
            dPath(ll:N:end, :), snapshotPtsT(ii), 1, ...
            'periodic', false);
    end
    c = deBoor(knotsS, nS, curveS, plotPtsS, 1, 'periodic', true);
    plot(c(:, 1), c(:, 2), 'k:', 'LineWidth', 1);
end

%% Finish plotting
hold off;