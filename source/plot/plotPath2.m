%% plotPath2
%
% Helper function to plot the opt path from the discretization and the true
% end point. 
%
% Input
%   dPath
%       The path
%   splineData
%       General information about the splines used.
%
function plotPath2(dPath,dend, splineData)

N = splineData.N;
Nt = splineData.Nt;
nT = splineData.nT;
knotsT = splineData.knotsT;

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

c0 = evalCurve(plotPtsS, d0, splineData);
c1 = evalCurve(plotPtsS, d1, splineData);
cend = evalCurve(plotPtsS, dend, splineData);
%% Setup plotting
hold on;
axis equal;
% axis tight;
% axis off;

%% Start  plotting
plot(c0(:, 1), c0(:, 2), 'k-', 'LineWidth', 1);
plot(c1(:, 1), c1(:, 2), 'k-', 'LineWidth', 1);
plot(cend(:, 1), cend(:, 2), 'r-', 'LineWidth', 2);


% Particle paths across time
for ii = 1:noLinesS
    for ll = Nt:-1:1
        lineT(ll, :) = evalCurve(linePtsS(ii), ...
                                 dPath((ll-1)*N+1:ll*N,:), splineData);
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
    c = evalCurve(plotPtsS, curveS, splineData);
    plot(c(:, 1), c(:, 2), 'k:', 'LineWidth', 1);
end

%% Finish plotting
hold off;