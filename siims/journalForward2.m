disp(mfilename);

%% Initial velocity experiment; wrap, propeller and corpus callosum
stepsTList = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60];
noStepsT = length(stepsTList);

%% Directory paths
% Uncomment for stand alone running
% dataDir = '~/diss/openprojects/h2_numerics/data/';
% prefixDir = '~/diss/openprojects/h2_numerics/journal/';
% plotDir = [prefixDir, 'final_plots/'];
workDir = [prefixDir, 'forward2/'];
shapesFile = [prefixDir, 'basic_shapes.mat'];
resultsFile = [workDir, 'results.mat'];

%% Load shapes
basicShapes = load(shapesFile);

splineDataList = basicShapes.splineDataList;
optionsList = basicShapes.optionsList;
d1List = basicShapes.d1List;
d2List = basicShapes.d2List;

% noCurves = length(splineDataList);
noCurves = 3;

%% Create directories
if ~exist(workDir, 'dir')
    mkdir(workDir);
end

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

%% Do calculations
velRelList = zeros(noStepsT, noCurves);

tic
for jj = 1:noCurves
    disp(['Curve ', num2str(jj)]);
    
    splineData = splineDataList{jj};
    splineData = setupQuadData(splineData);
    options = optionsList{jj};
    d1 = d1List{jj};
    d2 = d2List{jj};
    
    % Create initial velocity
    % Use geodesicBvp for that
    [optE, optPath] = geodesicBvp(d1, d2, splineData, 'options', options);
    v0 = pathVelocity(optPath, 0, splineData);
    % v0 = v0 ./ curveRiemH2Norm(d1, v0, splineData);
    
    for kk = noStepsT:-1:1
        disp(['stepsT=', num2str(stepsTList(kk))]);
        
        stepsT = stepsTList(kk);
        splineData.stepsT = stepsT;
        splineData.Nt = stepsT;
        splineData = constructKnots(splineData);
        splineData = setupQuadData(splineData);
        
        % Forward shooting
        dEnd = geodesicExp(d1, v0, splineData);
        
        filename = ['curve', num2str(jj), '_dend_stepsT', num2str(stepsT)]; 
        saveCurve(filename, dEnd, splineData, 'workDir', workDir);
        
        % Reconstruct initial velocity
        [optE, optPath] = geodesicBvp(d1, dEnd, splineData, ...
                                      'options', options);
        w0 = pathVelocity(optPath, 0, splineData);
        
        filename = ['curve', num2str(jj), '_w0_stepsT', num2str(stepsT)]; 
        saveCurve(filename, w0, splineData, 'workDir', workDir);
        
        velRelList(kk, jj) = curveRiemH2Norm(d1, v0-w0, splineData) / ...
                             curveRiemH2Norm(d1, v0, splineData);
    end
end
toc

save(resultsFile, 'velRelList');
disp('Calculations finished.');

%% Load results again
results = load(resultsFile);
velRelList = results.velRelList;

%% Plot parameters (in points)
lineWidth = 400;
figRelSize = 0.49;

%% Plot relative distances initial velocities
figRatio = 4/3;
sx = figRelSize * lineWidth;
sy = sx / figRatio;

handle = figure( 'PaperUnits', 'points', 'PaperSize', [sx, sy], ...
                 'Units', 'points', 'Position', [0, 0, sx, sy], ...
                 'Color', 'white' );
handle.Visible = 'off';

x = stepsTList;
y = velRelList;
loglog(x, y(:,1), 'k-x', 'MarkerSize', 4);
hold on;
semilogy(x, y(:,2), 'k--o', 'MarkerSize', 4);
semilogy(x, y(:,3), 'k:+', 'MarkerSize', 4);
hold off;
set(gca, 'FontSize', 6);
legend({'Circle to wrap', 'Propeller 3 to 4', 'Corpus callosum'}, ...
        'FontSize', 6);
    
figname = [ plotDir, 'geodesicForwardCompatible.eps' ];
export_fig(figname);