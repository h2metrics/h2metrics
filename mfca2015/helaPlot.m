prefix = prefixDir;
curvedir = [paramDir, 'spline/'];
savedir = [paramDir, 'finished_plots/'];

N = splineData.N;
nS = splineData.nS;
knotsS = splineData.knotsS;

%% Print numerical data for use in paper.
% Get list of curves
list_all = dir([prefix, curvedir, 'cell_*.mat']);
noCurves = length(list_all);

calibrationMatFile = matfile([prefix, calibrateFile]);
karcherInfoFile = matfile([prefix, paramDir, 'mean/mean_info.mat']);
pcaInfoFile = matfile([prefix, paramDir, 'pca/pca_results.mat']);

[~, splineData] = loadCurve([prefix, paramDir, 'mean/mean.mat']);

disp(' ');
disp('Summary of HeLa data for paper');
disp(['Number of curves: ', num2str(noCurves)]);
disp(['N=', num2str(splineData.N)]);
disp(['nS=', num2str(splineData.nS)]);
disp(['Typical energy=', num2str(typicalEnergy)]);
disp(['Average length=', num2str(calibrationMatFile.averageLength)]);
disp(['Coefficients in metric: ', num2str(splineData.a(1)), ' ', ...
    num2str(splineData.a(2)), ' ', num2str(splineData.a(3))]);
disp(['Weights for calibration: ', num2str(weightL2), ' ', ...
    num2str(weightH1), ' ', num2str(weightH2)]);

disp(' ');
disp('Karcher mean calculation');

karcherIterInfo = karcherInfoFile.info;
disp(['Number of iterations: ', num2str(length(karcherIterInfo))]);
    
karcherLastIter = karcherIterInfo(end);
disp(['Final objective function value: ', num2str(karcherLastIter.cost)]);
disp(['Final gradient norm value: ', num2str(karcherLastIter.gradnorm)]);

disp(' ');
disp('PCA from mean');

Lambda = pcaInfoFile.Lambda;
varExplained = cumsum(Lambda(1:6)') / sum(Lambda);
disp(['First 6 eigenvalues: ', num2str(Lambda(1:6)')]);
disp(['Explained variance: ', num2str(varExplained)]);

%% Directories setup
if ~exist([prefix, savedir], 'dir')
    mkdir([prefix, savedir]);
end

%% Plot 1, Some curves + karcher mean
loaddir = [paramDir, 'spline/'];
meanfile = [paramDir, 'mean/mean.mat'];

ind = [ 2, 5, 8, 13, 21, 38, 70, 71];
noCurves = length(ind);

% Get list of curves
list_all = dir([prefix, loaddir, 'cell_*.mat']);

% Load curves
dList = {};
for kk = noCurves:-1:1
    [d, ~] = loadCurve( list_all(ind(kk)).name, ...
                        'workdir', [prefix, loaddir] );
    dList{kk} = d;
end

% Load mean
[dMean, ~] = loadCurve( meanfile, 'workdir', prefix);

% Start plotting
handle = figure( 'Units', 'centimeters', 'Position', [0, 0, 12 4], ...
                 'Color', 'white' );
handle.Visible = 'off';
clf;
hold on;
axis equal;
axis off;

noPlotPtsS = 100;
plotPtsS = linspace(0, 2*pi, noPlotPtsS);

xmin = 0;
xmax = 0;
ymin = 0;
ymax = 0;

ax = [];

for jj = 2:-1:1
    for kk = 4:-1:1
        ax(4*(jj-1) + kk) = axes( 'Clipping', 'off', 'Position', ...
                                  [(kk-1)/6, (jj-1)/2, 1/6, 1/2] );
        c = deBoor( knotsS, nS, dList{4*(jj-1)+kk}, plotPtsS, 1, ...
                    'periodic', true );
        plot(c(:,1), c(:,2), 'k-', 'Clipping', 'off');
        
        xl = xlim();
        xmin = min(xmin, xl(1));
        xmax = max(xmax, xl(2));
        xlim([xmin, xmax]);
        
        yl = ylim();
        ymin = min(ymin, yl(1));
        ymax = max(ymax, yl(2));
        ylim([ymin, ymax]);
        
        axis equal;
        axis off;
        axis tight;
    end
end
xl = xlim();
xlim([xl(1)-0.05*(xl(2)-xl(1)), xl(2)+0.05*(xl(2)-xl(1))]);
yl = ylim();
ylim([yl(1)-0.05*(yl(2)-yl(1)), yl(2)+0.05*(yl(2)-yl(1))]);

linkaxes(ax);

% Now plot the mean(s)
axMean = axes( 'Clipping', 'off', ...
               'Position', [4/6, 0, 2/6, 1]);
c = deBoor( knotsS, nS, dMean, plotPtsS, 1, ...
            'periodic', true );
plot(c(:,1), c(:,2), 'k-', 'Clipping', 'off');
axis equal;
axis off;
axis tight;

xl = xlim();
xlim([xl(1)-0.05*(xl(2)-xl(1)), xl(2)+0.05*(xl(2)-xl(1))]);
yl = ylim();
ylim([yl(1)-0.05*(yl(2)-yl(1)), yl(2)+0.05*(yl(2)-yl(1))]);

% Finish up plotting
hold off;

figname = [ prefix, savedir, 'hela_1_cells.pdf' ];
export_fig(figname);

%% Plot 2, v1, principal directions 4x1
prefix = prefixDir;
loaddir = [paramDir, 'pca/'];
meanfile = [paramDir, 'mean/mean.mat'];

% Load mean
[dMean, ~] = loadCurve(meanfile, 'workdir', prefix);

dSpace = splineData.dSpace;
N = splineData.N;
nS = splineData.nS;
knotsS = splineData.knotsS;
quadData = setupQuadData(splineData);

% Get list of principal components
list_all = dir([prefix, loaddir, 'pc_*.mat']);
noPC = length(list_all);

% Load principal components
vList = {};
for jj = noPC:-1:1
    [v, ~] = loadCurve( list_all(jj).name, ...
                        'workdir', [prefix, loaddir] );
    vList{jj} = v;
end

% Setup plotting
handle = figure( 'Units', 'centimeters', 'Position', [0, 0, 32 8], ...
                 'Color', 'white' );
handle.Visible = 'off';
clf;
hold on;
axis equal;
axis tight;
axis off;

noPlotPtsS = 100;
noLinesS = 20;
plotPtsS = linspace(0, 2*pi, noPlotPtsS);
linePtsS = linspace(0, 2*pi, noLinesS+1);
linePtsS = linePtsS(1:end-1);

stepsT = splineData.stepsT;
noLambda = 3;

xmin = 0;
xmax = 0;
ymin = 0;
ymax = 0;

ax = [];

xp = [0, 1/4, 2/4, 3/4];
yp = [0, 0, 0, 0];
for jj = 4:-1:1
    v0 = vList{jj};
    
    disp(jj);
    
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
    
    ax(jj) = axes( 'Clipping', 'off', 'Position', ...
                   [xp(jj), yp(jj), 1/4, 1] );
    axis equal;
    axis off;
    hold on;

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
    
    % Setup limits
    xl = xlim();
    xmin = min(xmin, xl(1));
    xmax = max(xmax, xl(2));
    xlim([xmin, xmax]);
    
    yl = ylim();
    ymin = min(ymin, yl(1));
    ymax = max(ymax, yl(2));
    ylim([ymin, ymax]);

    hold off;
end

linkaxes(ax);

% Finish the plot
hold off;

figname = [ prefix, savedir, 'hela_2a_pc.pdf'];
export_fig(figname);

%% Plot 2, v2, principal directions - 3x2
prefix = prefixDir;
loaddir = [paramDir, 'pca/'];
meanfile = [paramDir, 'mean/mean.mat'];

% Load mean
[dMean, ~] = loadCurve(meanfile, 'workdir', prefix);

dSpace = splineData.dSpace;
N = splineData.N;
nS = splineData.nS;
knotsS = splineData.knotsS;
quadData = setupQuadData(splineData);

% Get list of principal components
list_all = dir([prefix, loaddir, 'pc_*.mat']);
noPC = length(list_all);

% Load principal components
vList = {};
for jj = noPC:-1:1
    [v, ~] = loadCurve( list_all(jj).name, ...
                        'workdir', [prefix, loaddir] );
    vList{jj} = v;
end

% Setup plotting
handle = figure( 'Units', 'centimeters', 'Position', [0, 0, 24 16], ...
                 'Color', 'white' );
handle.Visible = 'off';
clf;
hold on;
axis equal;
axis tight;
axis off;

noPlotPtsS = 100;
noLinesS = 20;
plotPtsS = linspace(0, 2*pi, noPlotPtsS);
linePtsS = linspace(0, 2*pi, noLinesS+1);
linePtsS = linePtsS(1:end-1);

stepsT = splineData.stepsT;
noLambda = 3;

xmin = 0;
xmax = 0;
ymin = 0;
ymax = 0;

ax = [];

xp = [0, 1/3, 2/3, 0, 1/3, 2/3];
yp = [1/2, 1/2, 1/2, 0, 0, 0];
for jj = 6:-1:1
    disp(jj);
    
    v0 = vList{jj};
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
    
    ax(jj) = axes( 'Clipping', 'off', 'Position', ...
                   [xp(jj), yp(jj), 1/3, 1/2] );
    hold on;

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
    
    % Setup limits
    xl = xlim();
    xmin = min(xmin, xl(1));
    xmax = max(xmax, xl(2));
    xlim([xmin, xmax]);
    
    yl = ylim();
    ymin = min(ymin, yl(1));
    ymax = max(ymax, yl(2));
    ylim([ymin, ymax]);

    axis equal;
    axis off;
    % axis tight;
    hold off;
end

% % Make limits slighly larger
% xl = xlim();
% xlim([xl(1)-0.05*(xl(2)-xl(1)), xl(2)+0.05*(xl(2)-xl(1))]);
% yl = ylim();
% ylim([yl(1)-0.05*(yl(2)-yl(1)), yl(2)+0.05*(yl(2)-yl(1))]);

linkaxes(ax);

% Finish the plot
hold off;

figname = [ prefix, savedir, 'hela_2b_pc.pdf'];
export_fig(figname);

%% Plot 3, random samples
loaddir = [paramDir, 'pca/'];

% Get random samples
list_all = dir([prefix, loaddir, 'random_sample_*.mat']);
noSamples = 6;

% Load samples
dList = {};
for jj = noSamples:-1:1
    [d, ~] = loadCurve( list_all(jj).name, ...
                        'workdir', [prefix, loaddir] );
    dList{jj} = d;
end

% Setup plotting
handle = figure( 'Units', 'centimeters', 'Position', [0, 0, 24 4], ...
                 'Color', 'white' );
handle.Visible = 'off';
clf;
hold on;
axis equal;
axis tight;
axis off;

xmin = 0;
xmax = 0;
ymin = 0;
ymax = 0;

ax = [];

for jj = noSamples:-1:1
    d0 = dList{jj};
    
    ax(jj) = axes( 'Clipping', 'off', 'Position', ...
                   [(jj-1)/6, 0, 1/6, 1] );
    axis equal;
    axis off;
    hold on;
    
    c0 = deBoor(knotsS, nS, d0, plotPtsS, 1, 'periodic', true);
    plot(c0(:, 1), c0(:, 2), 'k-', 'LineWidth', 1);
    
    % Setup limits
    xl = xlim();
    xmin = min(xmin, xl(1));
    xmax = max(xmax, xl(2));
    xlim([xmin, xmax]);
    
    yl = ylim();
    ymin = min(ymin, yl(1));
    ymax = max(ymax, yl(2));
    ylim([ymin, ymax]);
    
    hold off;
    
end

linkaxes(ax);

% Finish the plot
figname = [ prefix, savedir, 'hela_3_random_samples.pdf'];
export_fig(figname);

%% Plot 4, explained variance
% prefix = '~/diss/openprojects/h2_numerics/hela_cells/';
% loaddir = 'mean_aligned_csp_A3_B1_C6_n4_N12/';
% resultsfile = 'mean_alignment_results.mat';
% 
% pcaFile = matfile([prefix, loaddir, resultsfile]);
% 
% Lambda = pcaFile.Lambda;
% 
% x = 1:length(Lambda);
% y = cumsum(Lambda) / sum(Lambda);
% 
% handle = figure( 'Color', 'white', 'Units', 'centimeters', ...
%                  'Position', [0, 0, 20, 16] );
% plot(x, y, 'k-x', 'LineWidth', 1);
% grid on;
% 
% figname = [prefix, savedir, 'hela_4_explained_variance.pdf'];
% export_fig(figname);

%% Plot 5, 2d plot in first 2 principal components
prefix = prefixDir;
loaddir = [paramDir, 'spline/'];
resultsfile = [paramDir, 'pca/pca_results.mat'];
meanfile = [paramDir, 'mean/mean.mat'];

list_all = dir([prefix, loaddir, 'cell_*.mat']);
noCurves = length(list_all);

if ~isempty(maxNoCurves)
    noCurves = maxNoCurves;
end

% Load mean
[dMean, ~] = loadCurve( meanfile, 'workdir', prefix);
N = splineData.N;

% Load curves
dList = {};
for kk = noCurves:-1:1
    [d, ~] = loadCurve( list_all(kk).name, ...
                        'workdir', [prefix, loaddir] );
    dList{kk} = d;
end

pcaFile = matfile([prefix, resultsfile]);
pts2d = pcaFile.pts2d;

noPlotPtsS = 100;
plotPtsS = linspace(0, 2*pi, noPlotPtsS);

handle = figure( 'Color', 'white', 'Units', 'centimeters', ...
                 'Position', [0, 0, 20, 16] );
handle.Visible = 'off';
ax = axes('XTick', -3:2, 'YTick', -2:2);
hold on;

% Plot mean
dMean = curveCenter(dMean, splineData, quadData);
dAct = dMean ./ curveLength(dMean, splineData, quadData);

c = deBoor( knotsS, nS, dAct, plotPtsS, 1, 'periodic', true );
plot(c(:,1), c(:,2), 'b-', 'LineWidth', 2, 'Clipping', 'off');
plot(0, 0, 'b.');

% Plot 2d coordinates
plot(pts2d(1,:), pts2d(2,:), 'k.');

% Plot curve around each coordinate
options = struct( 'optTra', true, 'optRot', true, 'optShift', true, ...
                  'rigidUseComp', false, 'rigidGlobalRot', false );
for jj = noCurves:-1:1
    dTmp = rigidAlignment( {dMean, dList{jj}}, splineData, quadData, ...
                           'options', options );
    dAct = dTmp{2};
    dAct = curveCenter(dAct, splineData, quadData);
    dAct = dAct ./ curveLength(dAct, splineData, quadData);
    dAct = dAct + ones([N, 1]) * pts2d(:,jj)';
    
    c = deBoor(knotsS, nS, dAct, plotPtsS, 1, 'periodic', true);
    plot(c(:,1), c(:,2), 'k-', 'LineWidth', 0.5);
end

axis equal;
grid on;
box on;
hold off;

figname = [prefix, savedir, 'hela_5_2d_coordinates.pdf'];
export_fig(figname);


