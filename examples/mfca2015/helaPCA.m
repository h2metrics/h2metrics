%% Prerequisites
% splineData, quadData
% prefixDir, paramDir

disp('Do PCA analysis');

%% Parameters
noPC = 6;
noSamples = 10;
noLambda = 3;

N = splineData.N;
nS = splineData.nS;
knotsS = splineData.knotsS;
dSpace = splineData.dSpace;

%% Load curves
prefix = prefixDir;
loaddir = [paramDir, 'mean/'];
meanfile = [paramDir, 'mean/mean.mat'];
savedir = [paramDir, 'pca/'];
resultsfile = [prefix, savedir, 'pca_results.mat'];

list_all = dir([prefix, loaddir, 'path_*.mat']);
noCurves = length(list_all);

if ~isempty(maxNoCurves)
    noCurves = maxNoCurves;
end

% Load mean
[dMean, ~] = loadCurve( meanfile, 'workdir', prefix);

%% Directories setup
if ~exist([prefix, savedir], 'dir')
    mkdir([prefix, savedir]);
end

%% Load paths
dPathList = cell([noCurves, 1]);

for jj = noCurves:-1:1    
    filename = ['path_mean_to_', num2str(jj, '%02u'), '.mat'];
    [dPathTmp, ~, ~] = loadPathGamma( filename, ...
                                      'workdir', [prefix, loaddir] );
    dPathList{jj} = dPathTmp;
end

%% Do some PCA
G = metricMatrixH2(dMean, splineData, quadData);
rootG = sqrtm(G);

vList = {};
for jj = noCurves:-1:1
    vList{jj} = pathVelocity(dPathList{jj}, 0, splineData);
end

Sigma = zeros([N*dSpace, N*dSpace]);
for jj = noCurves:-1:1
    v = reshape(vList{jj}, [N*dSpace, 1]);
    Sigma = Sigma + v * v';
end

Sigma = 1./(noCurves-1) * rootG * Sigma * rootG;

[U, Lambda] = eig(Sigma);
Lambda = real(diag(Lambda));

%% Save principal components
for jj = 1:noPC
    v0 = rootG \ U(:,jj);
    v0 = sqrt(Lambda(jj)) * v0;
    v0 = reshape(v0, [N, dSpace]);
    
    filename = ['pc_', num2str(jj, '%02u.mat')];
    saveCurve( filename, v0, splineData, 'workdir', [prefix, savedir], ...
               'manualName', true );
end

%% Calculate 2d coordinates
V = U(:,1:2);
pts2d = zeros([2, noCurves]);
for jj = noCurves:-1:1
    v = reshape(vList{jj}, [N*dSpace, 1]);
    pts2d(:,jj) = V' * rootG * v;
    pts2d(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
end

%% Save all results
save( resultsfile, 'U', 'Lambda', 'Sigma', 'G', 'rootG', 'pts2d' );

%% Random sampling
disp('Calculate random samples');
noPlotPtsS = 100;
plotPtsS = linspace(0, 2*pi, noPlotPtsS);

for jj = noSamples:-1:1
    disp(jj);
    u = mvnrnd(zeros([N*dSpace,1]), diag(Lambda))';
    v0 = U * u;
    v0 = rootG \ v0;
    v0 = reshape(v0, [N, dSpace]);

    d1 = geodesicForward( dMean, dMean + v0 / splineData.stepsT, ...
                          splineData.stepsT, splineData, quadData, ...
                          'endpoint');
    d1 = curveCenter(d1, splineData, quadData);
    
    filename = ['random_sample_', num2str(jj, '%02u')];
    saveCurve( filename, d1, splineData, 'workdir', [prefix, savedir], ...
               'manualName', true );
           
    c1 = deBoor( knotsS, nS, d1, plotPtsS, 1, 'periodic', true);
    
    handle = figure;
    handle.Visible = 'off';
    
    plot(c1(:,1), c1(:,2), 'k-');
    
    figname = ['random_sample_', num2str(jj, '%02u.svg')];
    saveTightFigure(handle, [prefix, savedir, figname]);
end

%% Plot first noPC components
disp('Plot principal components');
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
stepsT = splineData.stepsT;

for jj = noPC:-1:1
    disp(jj);
    v0 = rootG \ U(:,jj);
    v0 = sqrt(Lambda(jj)) * v0;
    v0 = reshape(v0, [N, dSpace]);
    
    pcaPathP = geodesicForward( dMean, dMean + v0 / stepsT, ...
                                noLambda*stepsT+1, splineData, quadData);
    pcaPathM = geodesicForward( dMean, dMean - v0 / stepsT, ...
                                noLambda*stepsT+1, splineData, quadData);
    
    d0 = dMean;
    dSnapshots = zeros([N, dSpace, 2*noLambda]);
    for kk = noLambda:-1:1
        dSnapshots(:,:,kk) = pcaPathP(:,:,kk*stepsT+1);
        dSnapshots(:,:,noLambda+kk) = pcaPathM(:,:,kk*stepsT+1);
    end

    handle = figure;
    handle.Visible = 'off';
    clf;

    hold on;
    axis equal
    
    c0 = deBoor(knotsS, nS, d0, plotPtsS, 1, 'periodic', true);
    plot(c0(:, 1), c0(:, 2), 'k-', 'LineWidth', 2);

    % Particle paths across time
    for ii = 1:noLinesS
        for ll = size(pcaPathP, 3):-1:1
            linePT(ll, :) = deBoor(knotsS, nS, ...
                pcaPathP(:,:,ll), linePtsS(ii), 1, ...
                'periodic', true);
        end
        for ll = size(pcaPathM, 3):-1:1
            lineMT(ll, :) = deBoor(knotsS, nS, ...
                pcaPathM(:,:,ll), linePtsS(ii), 1, ...
                'periodic', true);
        end
        plot(linePT(:, 1), linePT(:, 2), 'k:', 'LineWidth', 1);
        plot(lineMT(:, 1), lineMT(:, 2), 'k:', 'LineWidth', 1);
    end

    % Snapshots of curve
    for ii = 1:size(dSnapshots, 3)
        c = deBoor( knotsS, nS, dSnapshots(:,:,ii), plotPtsS, 1, ...
                    'periodic', true);
        plot(c(:, 1), c(:, 2), 'k-', 'LineWidth', 1);
    end
    hold off;

    % Add some info
    someInfo = ['Principal direction', num2str(jj)];
    ylims = get(gca,'ylim');
    xlims = get(gca,'xlim');
    text(xlims(1),ylims(2)-0.1, someInfo);

    figname = [ 'principal_direction_', ...
                num2str(jj, '%02u'), '.svg'];
    saveTightFigure(handle, [prefix, savedir, figname]);
end

disp('PCA finished');
