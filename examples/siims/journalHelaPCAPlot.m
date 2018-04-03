disp(mfilename);

%% Plot principal directions
setList = {'hela_1d', 'hela_1p'};

for jj = 1:length(setList)
    setName = setList{jj};
    resultsFile = [prefixDir, 'hela/', setName, '.mat'];
    results = matfile(resultsFile);

    %% Extract things
    splineData = results.splineData;
    N = splineData.N;
    dSpace = splineData.dSpace;

    noCurves = results.noCurves;
    dMean = results.mean;
    meanPCList = results.meanPCList;

    %% Plot, principal directions 5x1
    noPC = 5;

    % Setup plotting
    handle = figure( 'Units', 'centimeters', 'Position', [0, 0, noPC*8 8], ...
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

    xp = (0:noPC-1) / noPC;
    yp = zeros(noPC, 1);

    for jj = noPC:-1:1
        disp(jj);  

        ax(jj) = axes( 'Clipping', 'off', 'Position', ...
                       [xp(jj), yp(jj), 1/noPC, 1] );
        axis equal;
        axis off;
        hold on;

        v = meanPCList{jj};
        plotGeodesic(dMean, v, splineData, [-3, 3], 7, 20);

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

    figname = [ plotDir, setName, '_', num2str(noPC), 'pc.eps'];
    export_fig(figname);
end

%% Plot mean
resultsFile = [prefixDir, 'hela/hela_1d.mat'];
results = matfile(resultsFile);
dMean = results.mean;
splineData = results.splineData;

% Setup plotting
lineWidth = 400;
figRelSize = 0.25;

figRatio = 3/3;
sx = figRelSize * lineWidth;
sy = sx / figRatio;
handle = figure( 'PaperUnits', 'points', 'PaperSize', [sx, sy], ...
                 'Units', 'points', 'Position', [0, 0, sx, sy], ...
                 'Color', 'white' );
handle.Visible = 'on';

axis equal;
axis tight;
axis off;

hold on;
noPlotPts = 300;
plotPts = linspace(0, 2*pi, noPlotPts);

c = evalCurve(plotPts, dMean, splineData);
plot(c(:,1), c(:,2), 'k-', 'LineWidth', 2);

figname = [ plotDir, 'hela_mean.eps'];
export_fig(figname);

%% 2D principal component plot
resultsFile = [prefixDir, 'hela/hela_1d.mat'];
results = matfile(resultsFile);

splineData = results.splineData;
N = splineData.N;
dSpace = splineData.dSpace;

noCurves = results.noCurves;
dList = results.dList;
dMean = results.mean;
pts2d = results.meanPts2d;
meanPCList = results.meanPCList;

% Now do plot
noPlotPtsS = 100;
plotPtsS = linspace(0, 2*pi, noPlotPtsS);

handle = figure( 'Color', 'white', 'Units', 'centimeters', ...
                 'Position', [0, 0, 20, 16] );
handle.Visible = 'on';
ax = axes('XTick', -2:3, 'YTick', -3:1);
hold on;

% Plot mean
dMean = curveCenter(dMean, splineData);
dAct = dMean ./ curveLength(dMean, splineData);

c = evalCurve(plotPtsS, dAct, splineData);
plot(c(:,1), c(:,2), 'b-', 'LineWidth', 2, 'Clipping', 'off');
plot(0, 0, 'b.');

% Plot 2d coordinates
plot(pts2d(1,:), pts2d(2,:), 'k.');

% Plot curve around each coordinate
options = struct( 'optTra', true, 'optRot', true, 'optShift', true, ...
                  'rigidUseComp', false, 'rigidGlobalRot', false );
for jj = noCurves:-1:1
    dTmp = rigidAlignment( {dMean, dList{jj}}, splineData, ...
                           'options', options );
    dAct = dTmp{2};
    dAct = curveCenter(dAct, splineData);
    dAct = dAct ./ curveLength(dAct, splineData);
    dAct = dAct + ones([N, 1]) * pts2d(:,jj)';
    
    c = evalCurve(plotPtsS, dAct, splineData);
    plot(c(:,1), c(:,2), 'k-', 'LineWidth', 0.5);
end

axis equal;
grid on;
box on;
hold off;

figname = [plotDir, 'hela_1d_pts2d.eps'];
export_fig(figname);