disp(mfilename);

%% Initial velocity experiment; corpus callosum data
NtList = [10, 15, 20, 30, 35, 40, 45, 50, 55, 60];
nTList = [1, 2, 3, 4];
noNt = length(NtList);
nonT = length(nTList);

%% Directory paths
% Uncomment for stand alone running
% dataDir = '~/diss/openprojects/h2_numerics/data/';
% prefixDir = '~/diss/openprojects/h2_numerics/journal/';
% plotDir = [prefixDir, 'final_plots/'];
workDir = [prefixDir, 'initvel/'];
shapesFile = [prefixDir, 'basic_shapes.mat'];
resultsFile = [workDir, 'results.mat'];
exampleInd = 3; % Corpus Callosum

%% Load shapes
basicShapes = load(shapesFile);

splineDataBase = basicShapes.splineDataList{exampleInd};
splineDataBase = setupQuadData(splineDataBase);
options = basicShapes.optionsList{exampleInd};
d1 = basicShapes.d1List{exampleInd};
d2 = basicShapes.d2List{exampleInd};

%% Create directories
if ~exist(workDir, 'dir')
    mkdir(workDir);
end

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

%% Do calculations
tic            
distList = zeros(noNt, nonT);
distEndList = zeros(noNt, nonT);

for ii = nonT:-1:1
    for jj = noNt:-1:1
        disp(['nT=', num2str(nTList(ii)), ' Nt=', num2str(NtList(jj))]);

        splineData = splineDataBase;
        splineData.Nt = NtList(jj);
        splineData.nT = nTList(ii);
        splineData.stepsT = NtList(jj);
        splineData = constructKnots(splineData);
        splineData = setupQuadData(splineData);

        [optE, optPath] = geodesicBvp( d1, d2, splineData, ...
                                       'options', options );
        distList(jj, ii) = sqrt(optE);

        fileName = 'path';
        savePath( fileName, optPath, splineData, 'workdir', workDir );

        stepsT = splineData.stepsT;
        v0 = pathVelocity(optPath, 0, splineData);
        dEnd = geodesicForward( d1, d1 + v0 ./ stepsT, stepsT, ...
                                splineData, 'endpoint' );

        optE2 = geodesicBvp(d2, dEnd, splineDataBase, 'options', options);
        distEndList(jj, ii) = sqrt(optE2);

        fileName = ['v0_nT', num2str(nTList(ii)), ...
                    '_Nt', num2str(NtList(jj))];
        saveCurve( fileName, v0, splineData, 'workdir', workDir);

        fileName = ['dEnd_nT', num2str(nTList(ii)), ...
                    '_Nt', num2str(NtList(jj))];
        saveCurve( fileName, dEnd, splineData, 'workdir', workDir);
    end
end
toc

save(resultsFile, 'distList', 'distEndList');
disp('Calculations finished.');

%% Load results again
results = matfile(resultsFile);
distList = results.distList;
distEndList = results.distEndList;

v0List = {};
dEndList = {};

for ii = length(nTList):-1:1
    for jj = length(NtList):-1:1
        fileName = ['v0_nT', num2str(nTList(ii)), ...
                    '_Nt', num2str(NtList(jj))];
        v0List{jj,ii} = loadCurve(fileName, 'workdir', workDir);
    
        fileName = ['dEnd_nT', num2str(nTList(ii)), ...
                    '_Nt', num2str(NtList(jj))];
        dEndList{jj,ii} = loadCurve(fileName, 'workdir', workDir);
    end
end

%% Plot parameters (in points)
lineWidth = 400;
figRelSize = 0.49;

%% Plot relative Riem. inner product difference for initial velocity
figRatio = 4/3;
sx = figRelSize * lineWidth;
sy = sx / figRatio;

handle = figure( 'PaperUnits', 'points', 'PaperSize', [sx, sy], ...
                 'Units', 'points', 'Position', [0, 0, sx, sy], ...
                 'Color', 'white' );
handle.Visible = 'off';

x = NtList(1:end-1);
y = zeros(noNt-1, nonT);
for ii = nonT:-1:1
    for jj = noNt-1:-1:1
        y(jj,ii) = curveRiemH2Norm(d1, v0List{jj,ii}-v0List{jj+1,ii}, ...
                                splineData) / ...
                curveRiemH2Norm(d1, v0List{jj+1,ii}, splineData);
    end
end
semilogy(x, y(:,1), 'k-x', 'MarkerSize', 4);
hold on;
semilogy(x, y(:,2), 'k--o', 'MarkerSize', 4);
semilogy(x, y(:,3), 'k:+', 'MarkerSize', 4);
semilogy(x, y(:,4), 'k-.*', 'MarkerSize', 4);
hold off;
set(gca, 'FontSize', 6);
legend({'$n_t=1$', '$n_t=2$', '$n_t=3$', '$n_t=4$'}, ...
        'FontSize', 6, 'Interpreter', 'latex');
    
figname = [ plotDir, 'initialVelocityRelative.eps' ];
export_fig(figname);

%% Plot geodesic distance of forward shooting to original endpoint
figRatio = 4/3;
sx = figRelSize * lineWidth;
sy = sx / figRatio;

handle = figure( 'PaperUnits', 'points', 'PaperSize', [sx, sy], ...
                 'Units', 'points', 'Position', [0, 0, sx, sy], ...
                 'Color', 'white' );
handle.Visible = 'off';

x = NtList;
y = distEndList ./ distList;
semilogy(x, y(:,1), 'k-x', 'MarkerSize', 4);
hold on;
semilogy(x, y(:,2), 'k--o', 'MarkerSize', 4);
semilogy(x, y(:,3), 'k:+', 'MarkerSize', 4);
semilogy(x, y(:,4), 'k-.*', 'MarkerSize', 4);
hold off;
set(gca, 'FontSize', 6);
legend({'$n_t=1$', '$n_t=2$', '$n_t=3$', '$n_t=4$'}, ...
        'FontSize', 6, 'Interpreter', 'latex');

figname = [ plotDir, 'initialVelocityDistEnd.eps' ];
export_fig(figname);
