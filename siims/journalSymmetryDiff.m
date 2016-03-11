disp(mfilename);

%% Symmetry of the geodesic distance
NList = [30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140];
NphiList = [10, 20, 30, 40];
indStart = [1, 2, 3, 4]; % Where to start N for each Nphi
nPhi = 3;

lineStyleList = {'k-x', 'k--o', 'k:+', 'k-.*'};
legendList = {'$N_\varphi=10$', '$N_\varphi=20$', ...
              '$N_\varphi=30$', '$N_\varphi=40$'};
noN = length(NList);% noN = 2;
noNphi = length(NphiList);% noNphi = 1;

%% Directory paths
% Uncomment for stand alone running
% dataDir = '~/diss/openprojects/h2_numerics/data/';
% prefixDir = '~/diss/openprojects/h2_numerics/journal/';
% plotDir = [prefixDir, 'final_plots/'];
workDir = [prefixDir, 'symm_diff/'];
shapesFile = [prefixDir, 'basic_shapes.mat'];
resultsFile = [workDir, 'results.mat'];

%% Load shapes
basicShapes = load(shapesFile);

splineDataList = basicShapes.splineDataList;
optionsList = basicShapes.optionsList;
d1List = basicShapes.d1List;
d2List = basicShapes.d2List;

noCurves = length(splineDataList);
%noCurves = 1;

%% Create directories
if ~exist(workDir, 'dir')
    mkdir(workDir);
end

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

%% Do calculations
distList = zeros(noN, noNphi, noCurves, 2);
distDiffList = zeros(noN, noNphi, noCurves);
for jj = 1:noCurves
    tic
    disp(['Curve ', num2str(jj)]);
    
    splineDataBase = splineDataList{jj};
    splineDataBase = setupQuadData(splineDataBase);
    options = optionsList{jj};
    options.optDiff = true;
    d1Base = d1List{jj};
    d2Base = d2List{jj};
    
    for kk = 1:noNphi
        disp([' Nphi=', num2str(NphiList(kk))]);
        
        for ll = indStart(kk):noN
            disp(['  N=', num2str(NList(ll))]);
            
            splineData = splineDataBase;
            splineData.N = NList(ll);
            splineData.Nphi = NphiList(kk);
            splineData.nPhi = nPhi;
            splineData = constructKnots(splineData);
            splineData = setupQuadData(splineData);
            
            d1 = curveSpline2Spline(d1Base, splineDataBase, splineData);
            d2 = curveSpline2Spline(d2Base, splineDataBase, splineData);
            
            % Computing the geodesic distance
            [optE1, optPath1, optGa1] = geodesicBvp(d1, d2, splineData, ...
                                            'options', options);
            [optE2, optPath2, optGa2] = geodesicBvp(d2, d1, splineData, ...
                                            'options', options);
                                        
            filename = ['curve', num2str(jj), '_1to2_Nphi', ...
                        num2str(NphiList(kk))];
            savePathGamma(filename, optPath1, optGa2, splineData, ...
                          'workDir', workDir);
            
            filename = ['curve', num2str(jj), '_2to1_Nphi', ...
                        num2str(NphiList(kk))];
            savePathGamma(filename, optPath2, optGa2, splineData, ...
                          'workDir', workDir);
            
            distList(ll, kk, jj, 1) = sqrt(optE1);
            distList(ll, kk, jj, 2) = sqrt(optE2);
            distDiffList(ll, kk, jj) = abs(sqrt(optE1)-sqrt(optE2)) / ...
                                       max(sqrt(optE1), sqrt(optE2));
        end
    end
    toc
end

save(resultsFile, 'distList', 'distDiffList');
disp('Calculations finished.');

%% Load results again
results = load(resultsFile);
distList = results.distList;
distDiffList = results.distDiffList;

%% Plot parameters (in points)
lineWidth = 400;
figRelSize = 0.49;

%% Plot relative distances between endpoints
figRatio = 4/3;
sx = figRelSize * lineWidth;
sy = sx / figRatio;

for jj = 1:noCurves
    handle = figure( 'PaperUnits', 'points', 'PaperSize', [sx, sy], ...
                     'Units', 'points', 'Position', [0, 0, sx, sy], ...
                     'Color', 'white' );
    handle.Visible = 'off';

    hold on;
    for kk = 1:noNphi
        x = NList(indStart(kk):noN);
        y = distDiffList(indStart(kk):noN, kk, jj);
        
        plot(x, y, lineStyleList{kk}, 'MarkerSize', 4);
    end
    hold off;
    set(gca, 'FontSize', 6);
    set(gca, 'YScale', 'log');
    legend(legendList(1:noNphi), 'Interpreter', 'latex', ...
           'FontSize', 6, 'Location', 'best');

    figname = [ plotDir, 'symmetricDistDiff', num2str(jj), '.eps' ];
    export_fig(figname);
end