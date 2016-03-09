%% Plot the Karcher mean from the Rohde paper
dataDir = '../data/';
sourceFile = [dataDir, 'source/hela_rohde/hela_rohde_mean.png'];
prefixDir = 'journal/';
plotDir = [prefixDir, 'final_plots/'];

% Load image file, use Otsu's method of thresholding, extract boundary
I = imread(sourceFile);
eff_luminance = graythresh(I);
BW = im2bw(I, eff_luminance);
[B, ~] = bwboundaries(BW, 4, 'noholes');
lengthComp = cellfun(@(t) size(t, 1), B);
[~, ind] = max(lengthComp);  
boundary = B{ind};
boundary = flip(boundary, 2); % Switch x- and y-coordinates

% Interpolate
auxN = 12;
auxnS = 4;
auxSplineData = constructEmptySplineData;
auxSplineData.N = auxN;
auxSplineData.nS = auxnS;
auxSplineData.dSpace = auxSplineData.dSpace;
auxSplineData = constructKnots(auxSplineData);
boundary = boundary(1:end-1,:); % Last point equals first
d0 = constructSplineApproximation(boundary, auxSplineData);

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

c = evalCurve(plotPts, d0, auxSplineData);
plot(c(:,1), c(:,2), 'k-', 'LineWidth', 2);

figname = [ plotDir, 'hela_rohde_mean.eps'];
export_fig(figname);


