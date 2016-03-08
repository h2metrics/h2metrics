loadDir = [ dataDir, 'source/hela_murphy/' ];
listAll = dir([loadDir, 'r*--1---2.dat.png']);

imInd = [1, 2, 3, 4];
noIm = 4;

%%
splineData = constructEmptySplineData;
splineData.N = 12;
splineData.nS = 4;
splineData.quadDegree = [6, 4];
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);

%%
dList = loadDataSet('hela_murphy', splineData, dataDir, 'ind', imInd);

%% Calculate plot dimensions
noPlotPts = 300;
plotPts = linspace(0, 2*pi, noPlotPts);

tmpfig = figure;

xdim = 0;
ydim = 0;
for jj = noIm:-1:1
    d0 = dList{jj};
    c0 = evalCurve(plotPts, d0, splineData);

    plot(c0(:,1), c0(:,2) , 'b', 'LineWidth', 2);
    
    xl = xlim();
    xdim = max(xdim, xl(2) - xl(1));
    yl = ylim();
    ydim = max(ydim, yl(2) - yl(1));
end
xdim = xdim / 2;
ydim = ydim / 2;
xdim = max(xdim, ydim);
ydim = max(xdim, ydim);

close(tmpfig);

%% Plot, some images
plotDir = [prefixDir, 'final_plots/'];

% Setup plotting
handle = figure( 'Units', 'centimeters', 'Position', [0, 0, noIm*8 8], ...
                 'Color', 'white' );
handle.Visible = 'off';
clf;
hold on;
axis equal;
axis tight;
axis off;

ax = [];

xp = (0:noIm-1) / noIm;
yp = zeros(noIm, 1);

for jj = noIm:-1:1
    disp(jj);  
    
    ax(jj) = axes( 'Clipping', 'off', 'Position', ...
                   [xp(jj), yp(jj), 1/noIm, 1] );
    
    sourceFile = [loadDir, listAll(jj).name];
    I = imread(sourceFile);
    
    d0 = dList{jj};
    [~, dCen] = curveCenter(d0, splineData);
    dCen = floor(dCen);
    c0 = evalCurve(plotPts, d0, splineData);

    plot(c0(:,1), c0(:,2) , 'b', 'LineWidth', 2);
    
    % Setup limits
    xlim([dCen(1) - xdim, dCen(1) + xdim]);
    ylim([dCen(2) - ydim, dCen(2) + ydim]);
    
    hold on;
    axis manual;
    imshow(I);
    plot(c0(:,1), c0(:,2) , 'b', 'LineWidth', 2);
    hold off;
   
end

% Finish the plot
hold off;

figname = [ plotDir, 'hela_images.eps'];
export_fig(figname);