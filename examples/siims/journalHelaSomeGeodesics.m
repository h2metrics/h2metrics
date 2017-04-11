%% Load curves
shapesFile = [prefixDir, 'basic_shapes.mat'];
basicShapes = load(shapesFile);
exampleInd = 4;

splineData = basicShapes.splineDataList{exampleInd};
splineData = setupQuadData(splineData);
options = basicShapes.optionsList{exampleInd};
d1 = basicShapes.d1List{exampleInd};
d2 = basicShapes.d2List{exampleInd};

splineDataLow = splineData;
splineDataLow.N = 12;
splineDataLow.nS = 4;
splineDataLow.Nphi = 6;
splineDataLow.nPhi = 4;
splineDataLow = constructKnots(splineDataLow);
splineDataLow = setupQuadData(splineDataLow);

d1Low = curveSpline2Spline(d1, splineData, splineDataLow);
d2Low = curveSpline2Spline(d2, splineData, splineDataLow);

%% Parametrized path N=12
options.optDiff = false;

[~, dPath1, gamma1, ~] = geodesicBvp(d1Low, d2Low, splineDataLow, ...
                                 'options', options);

%% Unparametrized path N=12
options.optDiff = true;

[~, dPath2, gamma2, ~] = geodesicBvp(d1Low, d2Low, splineDataLow, ...
                                 'options', options);
                             
%% Unparametrized path N=40
options.optDiff = true;

[~, dPath3, gamma3, ~] = geodesicBvp(d1, d2, splineData, ...
                                 'options', options);

%% Plot 1, some images
plotDir = [prefixDir, 'final_plots/'];
noIm = 3;

% Setup plotting
handle = figure( 'Units', 'centimeters', 'Position', [0, 0, noIm*4 4], ...
                 'Color', 'white' );
handle.Visible = 'off';
clf;
hold on;
axis equal;
axis tight;
axis off;

xmin = Inf;
xmax = -Inf;
ymin = Inf;
ymax = -Inf;

ax = [];

xp = (0:noIm-1) / noIm;
yp = zeros(noIm, 1);

noPlotPts = 300;
plotPts = linspace(0, 2*pi, noPlotPts);

dPathList = {dPath1, dPath2, dPath3};
dList = {d2Low, d2Low, d2};
gaList = {gamma1, gamma2, gamma3};
splineDataList = {splineDataLow, splineDataLow, splineData};
for jj = noIm:-1:1
    disp(jj);  
    
    ax(jj) = axes( 'Clipping', 'off', 'Position', ...
                   [xp(jj), yp(jj), 1/noIm, 1] );
               
    axis equal;
    axis off;
    hold on;
    
    plotPath(dPathList{jj}, splineDataList{jj});
    
    ga = gaList{jj};
    ga.phi = []; ga.alpha = [];
    d = curveApplyGamma(dList{jj}, ga, splineDataList{jj});
    if jj ~= 1
        plotCurve(d, splineDataList{jj}, 'b--');
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
   
end

linkaxes(ax);

% Finish the plot
hold off;

figname = [ plotDir, 'hela_geodesics.eps'];
export_fig(figname);