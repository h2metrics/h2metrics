disp(mfilename);

%% Directory paths
% Uncomment for stand alone running
% dataDir = '~/diss/openprojects/h2_numerics/data/';
% prefixDir = '~/diss/openprojects/h2_numerics/journal/';
% plotDir = [prefixDir, 'final_plots/'];
workDir = [prefixDir, 'noise/'];
shapesFile = [prefixDir, 'basic_shapes.mat'];
resultsFile = [workDir, 'results.mat'];

%% Load shapes - propeller
basicShapes = load(shapesFile);
exampleInd = 2;

splineData = basicShapes.splineDataList{exampleInd};
splineData = setupQuadData(splineData);
options = basicShapes.optionsList{exampleInd};
d1 = basicShapes.d1List{exampleInd};
d2 = basicShapes.d2List{exampleInd};

noise = [0, 2.^linspace(-13, -3, 11)];

%% Create directories
if ~exist(workDir, 'dir')
    mkdir(workDir);
end

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

%% H^2-continuity of the geodesic distance function
curveList = {'prop3', 'prop4'};
noCurves = 2;
noNoise = length(noise);

E_Riem = zeros(noCurves, noCurves, noNoise, noNoise, 4);
E_Flat = zeros(noCurves, noCurves, noNoise, noNoise, 4);

ind1List = [1, 1, 1];
ind2List = [1, 2, 2];
noise1List = [ones(1, noNoise); ones(1, noNoise); 1:noNoise];
noise2List = [1:noNoise; 1:noNoise; 1:noNoise];

tic
for jj = 1:length(ind1List)
    ind1 = ind1List(jj);
    ind2 = ind2List(jj);
    
    disp([curveList{ind1}, ' to ', curveList{ind2}]);

    for kk = 1:length(noise)
        noise1Ind = noise1List(jj,kk);
        noise2Ind = noise2List(jj,kk);
        
        noise1 = noise(noise1Ind);
        noise2 = noise(noise2Ind);

        disp(['  noise ', num2str(noise1), ...
              ' to noise ', num2str(noise2)]);

        d1 = loadDataSet('basic', splineData, '', ...
                     'curves', curveList{ind1}, 'noise', noise1);
        d2 = loadDataSet('basic', splineData, '', ...
                     'curves', curveList{ind2}, 'noise', noise2);
        
        [optE, optPath, optGa, info] = geodesicBvp(d1, d2, ...
            splineData, 'options', options);

        fileName = [curveList{ind1}, '_n', num2str(noise1), ...
            '_2_', curveList{ind2}, '_n', num2str(noise2) ];
        savePath(fileName, optPath, splineData, 'workdir', workDir);

        E_Riem(ind1, ind2, noise1Ind, noise2Ind,4) = optE;
        [~, comp] = pathRiemH2Energy(optPath, ...
            splineData, 'a', [1 1 1]);
        E_Riem(ind1, ind2, noise1Ind, noise2Ind,1:3) = comp;
        [E_Flat(ind1, ind2, noise1Ind, noise2Ind,4), comp] = ...
            curveFlatH2Norm(d1-d2, splineData, 'a', [1 1 1]);
        E_Flat(ind1, ind2, noise1Ind, noise2Ind,1:3) = comp;
    end
end
toc

save(resultsFile, 'E_Riem', 'E_Flat');
disp('Calculations finished.');

%% Load results again
resultsFile = load(resultsFile);
E_Flat = resultsFile.E_Flat;
E_Riem = resultsFile.E_Riem;
D_Riem = sqrt(E_Riem);

%% Plot parameters (in points)
lineWidth = 400;
figRelSize = 0.49;

%% Plot
figRatio = 4/3;
sx = figRelSize * lineWidth;
sy = sx / figRatio;
handle = figure( 'PaperUnits', 'points', 'PaperSize', [sx, sy], ...
                 'Units', 'points', 'Position', [0, 0, sx, sy], ...
                 'Color', 'white' );
handle.Visible = 'off';

y1 = D_Riem(1,1,1,:,4);
y2 = abs(D_Riem(1,2,1,:,4) - D_Riem(1,2,1,1,4)) / D_Riem(1,2,1,1,4);
dist = diag(squeeze(D_Riem(1,2,:,:,4)));
y3 = abs(dist - D_Riem(1,2,1,1,4)) / D_Riem(1,2,1,1,4);

loglog(noise, squeeze(y1), 'k-o', 'MarkerSize', 4);
hold on
loglog(noise, squeeze(y2), 'k--x', 'MarkerSize', 4);
loglog(noise, squeeze(y3), 'k-.d', 'MarkerSize', 4);
hold off
set(gca, 'XDir', 'reverse')
set(gca, 'FontSize', 6);
legend({'$c_0$ to $c_0+\varepsilon.n$', ...
        '$c_0$ to $c_1+\varepsilon.n$', ...
        '$c_0+\varepsilon.n$ to $c_1 +\varepsilon.n$'}, ...
        'Interpreter', 'latex', 'FontSize', 6, 'Location', 'southwest');
legend();

figname = [ plotDir, 'noisePropeller.eps' ];
export_fig(figname);