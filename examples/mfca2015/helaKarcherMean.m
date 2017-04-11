%% Prerequisites
% splineData, quadData, quadDataTensor
% prefixDir, paramDir

%% Parameters
options = struct( 'optDiff', false, ...
                  'optTra', true, ...
                  'optRot', true, ...
                  'optShift', true, ...
                  'tolFun', 1e-3, ...
                  'tolX', 1e-4, ...
                  'display', 'off', ... % 'iter-detailed'
                  'maxIter', 300, ...
                  'rigidUseComp', true, ...
                  'rigidGlobalRot', true, ...
                  'karcherTolGradNorm', 1e-3, ...
                  'karcherMaxIter', 50 );

%% Load splineData
prefix = prefixDir;
loaddir = [paramDir, 'spline/'];
savedir = [paramDir, 'mean/'];
meanguess = meanGuessFile;

list_all = dir([prefix, loaddir, 'cell_*.mat']);
noCurves = length(list_all);

if ~isempty(maxNoCurves)
    noCurves = maxNoCurves;
end

%% Directory setup
if ~exist([prefix, savedir], 'dir')
    mkdir([prefix, savedir]);
end

%% Load all curves
dList = {};

for kk = noCurves:-1:1
    [d, ~] = loadCurve(list_all(kk).name, 'workdir', [prefix, loaddir]);
    dList{kk} = d;
end

if ~isempty(meanguess) && exist([prefix, meanguess], 'file')
    [meanInit, ~] = loadCurve(meanguess, 'workdir', prefix);
else
    meanInit = [];
end

%% Compute Karcher mean
disp('Computing the Karcher mean');
[dMean, info] = karcherMeanManopt( dList, ...
    splineData, quadData, quadDataTensor, 'options', options, ...
    'meanInit', meanInit );
disp('Karcher mean computed');

dPathList = info(end).dPathList;
gaList = info(end).gaList;

%% Save the Karcher mean
filename = 'mean.mat';
saveCurve( filename, dMean, splineData, 'workdir', [prefix, savedir], ...
           'manualName', true );

filename = [prefix, savedir, 'mean_info.mat'];
save(filename, 'info', 'options');

%% Save all paths
for kk = noCurves:-1:1
    dPath = dPathList{kk};
    ga = gaList{kk};
    
    filename = ['path_mean_to_', num2str(kk, '%02u'), '.mat'];
    savePathGamma( filename, dPath, ga, splineData, ...
                   'workdir', [prefix, savedir], 'manualName', true );
end

%% Plot the geodesic paths
disp('Plotting geodesic paths');

N = splineData.N;
nS = splineData.nS;
knotsS = splineData.knotsS;

noPlotPtsS = 100;
noPlotPtsT = 100;
noSnapshotsT = 5;
noLinesS = 20;

plotPtsS = linspace(0, 2*pi, noPlotPtsS+1);
plotPtsT = linspace(0, 1, noPlotPtsT);
linePtsS = linspace(0, 2*pi, noLinesS+1);
linePtsS = linePtsS(1:end-1);
snapshotPtsT = linspace(0, 1, noSnapshotsT + 2);
snapshotPtsT = snapshotPtsT(2:end-1);
        
for jj = noCurves:-1:1
    dPathAct = dPathList{jj};
    d0 = dPathAct(1:N, :);
    d1 = dPathAct(end-N+1:end, :);

    c0 = deBoor(knotsS, nS, d0, plotPtsS, 1, 'periodic', true);
    c1 = deBoor(knotsS, nS, d1, plotPtsS, 1, 'periodic', true);

    handle = figure;
    handle.Visible = 'off';
    clf;

    hold on;
    axis equal
    plot(c0(:, 1), c0(:, 2), 'k-', 'LineWidth', 1);
    plot(c1(:, 1), c1(:, 2), 'k-', 'LineWidth', 1);

    % Particle paths across time
    for ii = 1:noLinesS
        for ll = splineData.Nt:-1:1
            lineT(ll, :) = deBoor(knotsS, nS, ...
                dPathAct((ll-1)*N+1:ll*N,:), linePtsS(ii), 1, ...
                'periodic', true);
        end
        b = deBoor(splineData.knotsT, splineData.nT, lineT, ...
            plotPtsT, 1, 'periodic', false);
        plot(b(:, 1), b(:, 2), 'k-', 'LineWidth', 1);
    end

    % Snapshots of curve
    curveS = [];
    for ii = 1:noSnapshotsT
        for ll = N:-1:1
            curveS(ll, :) = deBoor(splineData.knotsT, splineData.nT, ...
                dPathAct(ll:N:end, :), snapshotPtsT(ii), 1, ...
                'periodic', false);
        end
        c = deBoor(knotsS, nS, curveS, plotPtsS, 1, 'periodic', true);
        plot(c(:, 1), c(:, 2), 'k:', 'LineWidth', 1);
    end

    % Add some info
    someInfo = ['Mean to ', num2str(jj), ...
                ', distance ', num2str(info(end).enList(jj))];
    ylims = get(gca,'ylim');
    xlims = get(gca,'xlim');
    text(xlims(1),ylims(2)-0.1, someInfo);

    figname = [ 'path_mean_to_', ...
                num2str(jj, '%02u'), '.svg'];
    saveTightFigure(handle, [prefix, savedir, figname]);
end
