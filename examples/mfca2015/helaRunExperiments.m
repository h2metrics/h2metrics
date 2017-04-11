%% Setup parameters
splineData = constructEmptySplineData;
splineData.dSpace = 2;
splineData.N = 12;
splineData.nS = 4;
splineData.Nt = 10 + 2;
splineData.nT = 2;
splineData.Nphi = 5;
splineData.nPhi = 3;
splineData.quadDegree = [8 4];
splineData.noInterpolS = 4 * splineData.N;
splineData.a = [];
splineData.stepsT = 50;
splineData = constructKnots(splineData);
[quadData, quadDataTensor] = setupQuadData(splineData);

typicalEnergy = 100;
weightL2 = 3; % They should sum to 1
weightH1 = 1;
weightH2 = 6    ;

maxNoCurves = [];

%% Directory paths
prefixDir = '~/diss/openprojects/h2_numerics/hela_cells/';
sourceDir = 'murphy_source_images/';
splineDir = 'spline_n4_N12/';
splineCspDir = 'spline_csp_n4_N12/';
calibrateFile = 'calibration_n4_N12.mat';
paramDir = 'A3_B1_C6_n4_N12/';
meanGuessFile = [paramDir, 'spline/cell_24_n4_N12.mat'];

%% Run things
% helaFindCurves;

%% 
% helaCalibrateMetric;

%% 
helaApplyCalibration;

%%
helaKarcherMean;

%%
helaPCA;

%%
helaPlot;