%% Prerequisites
% splineData
% prefixDir, sourceDir, splineDir, splineCspDir

%% Setup interpolation data
% This spline is used to interpolate the boundary
auxN = 12;
auxnS = 4;
auxSplineData = constructEmptySplineData;
auxSplineData.N = auxN;
auxSplineData.nS = auxnS;
auxSplineData.dSpace = splineData.dSpace;
auxSplineData = constructKnots(auxSplineData);

% This spline data has interpolation parametrs set; don't use for saving
splineDataInterpol = splineData;
splineDataInterpol.noInterpolS = 12 * splineData.N;
splineDataInterpol = constructKnots(splineDataInterpol);
quadDataInterpol = setupQuadData(splineDataInterpol);

%% Prepare list of filenames
prefix = prefixDir;
loaddir = sourceDir;
savedir = splineDir;
cspsavedir = splineCspDir;

list_all = dir([prefix, loaddir, '*.png']);
noCurves = length(list_all);

if ~isempty(maxNoCurves)
    noCurves = maxNoCurves;
end

list = cell([noCurves, 1]);
for kk = 1:noCurves
    list{kk} = list_all(kk).name;
end

%% Check if directories exist
if ~exist([prefix, savedir], 'dir')
    mkdir([prefix, savedir]);
end

if ~exist([prefix, cspsavedir], 'dir')
    mkdir([prefix, cspsavedir]);
end

%% Do main work
disp('Start boundary extraction.');
for jj = 1:noCurves
    % Use Otsu's method of thresholding
    I = imread([prefix, loaddir, list{jj}]);
    eff_luminance = graythresh(I);
    disp([jj, eff_luminance]);

    BW = im2bw(I, eff_luminance);
    % BW = bwareaopen(BW, 12, 4);
    [B, ~] = bwboundaries(BW, 4, 'noholes');
    lengthComp = cellfun(@(t) size(t, 1), B);
    [~, ind] = max(lengthComp);  
      
    % Extract boundary
    boundary = B{ind};
    boundary = flip(boundary, 2); % Switch x- and y-coordinates
    
    % Do interpolation; interpolate first with auxSplineData, then
    % upsample to higher order
    boundary = boundary(1:end-1,:); % Last point equals first
    d0 = constructSplineApproximation(boundary, auxSplineData);
    dPts = deBoor(auxSplineData.knotsS, auxSplineData.nS, d0, ...
        splineDataInterpol.interpolS, 1, 'periodic', true);
    d0 = quadDataInterpol.B_interpolS \ dPts;
    
    % Save curve
    cellname = ['cell_', num2str(jj, '%02u')];
    filename = cellname;
    saveCurve(filename, d0, splineData, 'workdir', [prefix, savedir]);
    
    % Reparametrize to constant speed
    dCsp = curveReparamConstSpeed(d0, ...
        splineDataInterpol, quadDataInterpol);
    
    % Save curve
    filename = cellname;
    saveCurve(filename, dCsp, splineData, 'workdir', [prefix, cspsavedir]);
    
    
    % Superimpose boundary over cell -- curve
    handle = figure(2);
    handle.Visible = 'off';
    plot(boundary(:,1), boundary(:,2) , 'b', 'LineWidth', 2);
    hold on;
    axis manual;
    imshow(I);
    plot(boundary(:,1), boundary(:,2) , 'b', 'LineWidth', 2);
    hold off;
    
    % Superimpose interpolation over leaf
    noPlotPoints = size(boundary, 1);
    plotPoints = linspace(0, 2*pi, noPlotPoints);
    boundary2 = deBoor( splineData.knotsS, splineData.nS, d0, ...
                        plotPoints, 1, 'periodic', true);
    hold on;
    plot(boundary2(:,1), boundary2(:,2) , 'r', 'LineWidth', 2);
    hold off;
    
    % Save figure
    figname = [prefix, savedir, cellname, '.jpg'];
    saveTightFigure(handle, figname);
    
   
    % Superimpose boundary over leaf -- constant speed curve
    handle = figure(3);
    handle.Visible = 'off';
    plot(boundary(:,1), boundary(:,2) , 'b', 'LineWidth', 2);
    hold on;
    axis manual;
    imshow(I);
    plot(boundary(:,1), boundary(:,2) , 'b', 'LineWidth', 2);
    hold off;
    
    % Superimpose const speed interpolation over leaf
    noPlotPoints = size(boundary, 1);
    plotPoints = linspace(0, 2*pi, noPlotPoints);
    boundary3 = deBoor( splineData.knotsS, splineData.nS, dCsp, ...
                        plotPoints, 1, 'periodic', true);
    hold on;
    plot(boundary3(:,1), boundary3(:,2) , 'r', 'LineWidth', 2);
    hold off;
    
    % Save figure
    figname = [prefix, cspsavedir, cellname, '.jpg'];
    saveTightFigure(handle, figname);
end
disp('Boundary extraction finished.');
