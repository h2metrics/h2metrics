%% 
% Create bumps for paper

%% Set paths
clear;
addpath(genpath('./source'));
addpath('./lib/varifolds');
addpath('./lib/hanso2_2');
addpath('./lib/export_fig');
dataDir = '../../data';

resultsFile = 'paperBumps.mat';
expName = 'bumps3';

set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters');

%% Setup splineData
splineData = constructSplineData;
splineData.curveClosed = 0;
splineData.N = 120;
splineData.nS = 3;
splineData.Nt = 10;
splineData.nT = 2;

splineData.quadDegree = [6, 4];

varData = constructEmptyVarData(splineData);
varData.kernelGrass = 'gaussian_oriented';
varData.kernelSizeGeom = 0.1;
varData.kernelSizeGrass = 0.3;
varData.noPts = 3 * splineData.N;
varData.pts = [];
splineData.varData = varData;

splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);

splineData.scaleInv = 0;

aMap = containers.Map();
aMap('a_0_0_1e-3_10_1e-1') = [0 0 1e-3 10 1e-1];
aMap('a_0_0_1e-3_10_1') = [0 0 1e-3 10 1];
aMap('a_0_0_1e-3_10_1e1') = [0 0 1e-3 10 1e1];
aMap('a_0_0_1e-3_10_1e2') = [0 0 1e-3 10 1e2];
aMap('a_0_0_1e-3_10_1e-2') = [0 0 1e-3 10 1e-2];
aMap('a_0_0_1e-3_10_1e-3') = [0 0 1e-3 10 1e-3];
aMap('a_0_0_1_0_0') = [0 0 1 0 0];
aMap('a_1_0_1e-3_0_0') = [1 0 1e-3 0 0];

nameList = keys(aMap);
nameList = {'a_0_0_1e-3_10_1e-2', ...
            'a_0_0_1e-3_10_1e-1', ...
            'a_0_0_1e-3_10_1', ...
            'a_0_0_1e-3_10_1e1'};
nameList = {'a_0_0_1e-3_10_1e-3'};
nameList = {'a_0_0_1_0_0'};
nameList = {'a_1_0_1e-3_0_0'};
nameList = {'a_0_0_1e-3_10_1'};
nameList = keys(aMap);

%% Create the bumps
d1 = createBump(pi-5*pi/10, 4*pi/5, 4, splineData);
d2 = createBump(pi+5*pi/10, 4*pi/5, 4, splineData);

d1 = curveReparamConstSpeed(d1, splineData);
d2 = curveReparamConstSpeed(d2, splineData);

d1 = d1 / curveLength(d1, splineData) * 2*pi;
d2 = d2 / curveLength(d2, splineData) * 2*pi;

dInit = linearPath(d1, d2, splineData);

%% Create optimization options
options = struct( 'optDiff', true, ...
                  'optTra', false, ...
                  'optRot', false, ...
                  'useVarifold', true, ...
                  'varLambda', [], ...
                  'useMultigrid', true, ...
                  'useAugmentedLagrangian', false, ...
                  'hansoNvec', 500, ...
                  'hansoMaxIt', 1000, ...
                  'hansoCpuMax', 600, ...
                  'hansoNormTol', 1e-2, ...
                  'hansoPrtLevel', 2 );
mgOptions = options;
mgOptions.useMultigrid = false;
mgOptions.hansoMaxIt = 1000;
mgOptions.hansoNormTol = 1e-1;
options.mgOptions = mgOptions;

optQuick = options;
optQuick.varLambda = 100;
optQuick.useMultigrid = false;
optQuick.hansoMaxIt = 50;
optQuick.hansoPrtLevel = 0;

sdQuick = splineData;
sdQuick.N = splineData.N / 2;
sdQuick.varData.noPts = 3*sdQuick.N;
sdQuick.varData.pts = [];
sdQuick = constructKnots(sdQuick);
sdQuick = setupQuadData(sdQuick);

%% Load previous results
d1Map = containers.Map();
d2Map = containers.Map();
dInitMap = containers.Map();
splineDataMap = containers.Map();
optionsMap = containers.Map();
optPathMap = containers.Map();
optGaMap = containers.Map();
infoMap = containers.Map();

% Data for a geodesic calculation
% d1, d2, dInit
% splineData, options
% optPath, optGa, info

if exist(resultsFile, 'file')
    load(resultsFile, 'd1Map', 'd2Map', 'dInitMap', ...
         'splineDataMap', 'optionsMap', 'optPathMap', ...
         'optGaMap', 'infoMap');
end

%% Do the computations
for jj = 1:length(nameList)
    % Set right parameters
    splineData.a = aMap(nameList{jj});
    sdQuick.a = aMap(nameList{jj});
    
    d1q = curveSpline2Spline(d1, splineData, sdQuick);
    d2q = curveSpline2Spline(d2, splineData, sdQuick);
    dInitQuick = pathSpline2Spline(dInit, splineData, sdQuick);
    %dInitq = linearPath(d1q, d2q, sdQuick);
    
    varEnGoal = 1e-3; % My target for the squared varifold distance
    [distEnEst, ~, ~, ~] = geodesicBvp(d1q, d2q, sdQuick, ...
                                       optQuick, 'initPath', dInitQuick);
    
    % My estimate for lambda
    lambda = sqrt(distEnEst / varEnGoal);
    disp([num2str(distEnEst), ' ', num2str(sqrt(lambda))]);
    
    options.varLambda = lambda;
    options.mgOptions.varLambda = lambda;
    
    [~, optPath, optGa, info] = geodesicBvp(d1, d2, splineData, ...
                                         options, 'initPath', dInit);

    d1Map(nameList{jj}) = d1;
    d2Map(nameList{jj}) = d2;
    dInitMap(nameList{jj}) = dInit;
    splineDataMap(nameList{jj}) = splineData;
    optionsMap(nameList{jj}) = options;
    optPathMap(nameList{jj}) = optPath;
    optGaMap(nameList{jj}) = optGa;
    infoMap(nameList{jj}) = info;

    analyzeBvpResults(d1, d2, optPath, splineData, options, info);
                                      
    %videoPath(['bump_', num2str(jj), '.avi'], optPath, ...
    %          splineData, 'lineStyle', 'colour', 'lineWidth', 2);
end

%% Save data
save(resultsFile, 'd1Map', 'd2Map', 'dInitMap', 'splineDataMap', ...
     'optionsMap', 'optPathMap', 'optGaMap', 'infoMap');

%% Load data
clearvars -except 'dataDir' 'resultsFile' 'expName' 'nameList'
load(resultsFile, 'd1Map', 'd2Map', 'dInitMap', 'splineDataMap', ...
     'optionsMap', 'optPathMap', 'optGaMap', 'infoMap');

%% Plots
for jj = 1:length(nameList)
    d1 = d1Map(nameList{jj});
    d2 = d2Map(nameList{jj});
    dInit = dInitMap(nameList{jj});
    splineData = splineDataMap(nameList{jj});
    options = optionsMap(nameList{jj});
    optPath = optPathMap(nameList{jj});
    optGa = optGaMap(nameList{jj});
    info = infoMap(nameList{jj});
    
    analyzeBvpResults(d1, d2, optPath, splineData, options, info);
    
    plotPathRow(optPath, splineData, 'noPtsT', 5, ...
        'lineStyle', 'colour', 'lineWidth', 2);
    
    export_fig([expName, '_', nameList{jj}, '.pdf']);
    
    videoPath([expName, '_', nameList{jj}, '.avi'], optPath, ...
              splineData, 'lineStyle', 'colour', 'lineWidth', 2);
    
    close;
end