%%
% Assume tRange = [tmin, tmax] with tmin <= 0, tmax >= 0
% noSnapshots should be odd
% (noSnapshots-1) / 2 lie on ne
function plotGeodesic(d, v, splineData, tRange, noSnapshots, noParticles)

dSpace = splineData.dSpace;
N = splineData.N;
stepsT = splineData.stepsT;

%% Calculate snapshots of curve
tmin = tRange(1);
tmax = tRange(2);
if tmin < 0 && tmax > 0
    noSnapP = round((noSnapshots-1) / 2);
    noSnapM = round((noSnapshots-1) / 2);
elseif tmin < 0
    noSnapP = 0;
    noSnapM = noSnapshots - 1;
elseif tmax > 0
    noSnapP = noSnapshots - 1;
    noSnapM = 0;
else
    noSnapP = 0;
    noSnapM = 0;
end

dSnapshots = zeros([N, dSpace, noSnapP+noSnapM+1]);
dPath = zeros([N, dSpace, (noSnapP+noSnapM)*stepsT + 1]);

if noSnapM > 0
    dPathM = geodesicForward( d, d - v / stepsT, noSnapM * stepsT + 1, ...
                              splineData);
	dPath(:,:,1:noSnapM*stepsT) = flip(dPathM(:,:,2:end), 3);
    dSnapshots(:,:,1:noSnapM) = dPathM(:,:,stepsT+1:stepsT:end);
end

dPath(:,:,noSnapM*stepsT+1) = d;
dSnapshots(:,:,noSnapM+1) = d;

if noSnapP > 0
    dPathP = geodesicForward( d, d + v / stepsT, noSnapP * stepsT + 1, ...
                              splineData);
    dPath(:,:,noSnapM*stepsT+2:end) = dPathP(:,:,2:end);
    dSnapshots(:,:,noSnapM+2:end) = dPathP(:,:,stepsT+1:stepsT:end);
end

hold on;

%% Plot snapshots
noPlotPtsS = 100;
plotPtsS = linspace(0, 2*pi, noPlotPtsS);

% Snapshots of curve
for jj = 1:size(dSnapshots, 3)
    c = evalCurve(plotPtsS, dSnapshots(:,:,jj), splineData);
    if jj <= noSnapM
        plot(c(:, 1), c(:, 2), 'k--', 'LineWidth', 1);
    else
        plot(c(:, 1), c(:, 2), 'k-', 'LineWidth', 1);
    end
end

%% Plot curve itself
c = evalCurve(plotPtsS, d, splineData);
plot(c(:, 1), c(:, 2), 'k-', 'LineWidth', 2);

%% Plot lines
noLinesS = noParticles;
linePtsS = linspace(0, 2*pi, noLinesS+1);
linePtsS = linePtsS(1:end-1);

% Particle paths across time
for jj = 1:noLinesS
    for ll = size(dPath, 3):-1:1
        linePt(ll, :) = evalCurve(linePtsS(jj), dPath(:,:,ll), splineData);
    end
    plot(linePt(:, 1), linePt(:, 2), 'k:', 'LineWidth', 1);
end

hold off;
axis equal;
