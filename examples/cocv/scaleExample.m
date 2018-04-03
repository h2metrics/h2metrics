%% 
%Scale

%% Precaution
clear all;

%% Set paths
% This should be done in a setup script
addpath('./source');
addpath('./source/datasets');
addpath('./source/save');
addpath('./source/metric');
addpath('./source/plot');
addpath('./source/riemgeom');
addpath('./source/setup');
addpath('./source/utility');
addpath('./examples');
addpath('./lib/varifolds');
addpath('./lib/hanso2_2');
addpath('./lib/export_fig');
dataDir = '../../data';

%% Semantic parameters
charLength = 2*pi;  % Length of an average curve
relKernelSize = 2;  % How many knots to each side should the kernel see?
varEnGoal = 1e-3;   % Varifold norm squared should be this much smaller
                    % than the geodesic distance squared
a = [0 1 5e-4 0 0];

%% Setup splineData
splineData = constructSplineData;
splineData.curveClosed = 1;
splineData.N = 120;
splineData.nS = 3;
splineData.Nt = 10;
splineData.nT = 2;
splineData.scaleInv = 1;

splineData.quadDegree = [6, 4];

varData = constructEmptyVarData(splineData);
varData.kernelGrass = 'gaussian_unoriented';
varData.kernelSizeGeom = []; % charLength / splineData.N * relKernelSize;
varData.kernelSizeGrass = 0.3;
varData.noPts = 3 * splineData.N;
varData.pts = [];
splineData.varData = varData;

splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);

if ~splineData.scaleInv
    toScale = [(2*pi/charLength)^3, 2*pi/charLength, ...
               (2*pi/charLength)^(-1), 2*pi/charLength, 2*pi/charLength];
    splineData.a = a .* toScale;
else
    toScale = [(2*pi)^3, 2*pi, (2*pi)^(-1), 2*pi, 2*pi];
    splineData.a = a .* toScale;
end

%% Load some curves
dList = loadDataSet('surrey_fish', splineData, dataDir,...
                     'constspeed', 'ind', 1:20);

%d_list = loadDataSet('corpus_callosum_tilak', splineData, dataDir,...
%    'class','ad','noPlot','constspeed','ind',1:10);

noCurves = length(dList);
dListUnit = cell(noCurves, 1);
for ii = 1:noCurves
    ell = curveLength(dList{ii}, splineData);
    dListUnit{ii} = dList{ii} / ell * 2*pi;
end

%% Choose example class
% Corpus Callosi
d1orig = dListUnit{2};
d2orig = dListUnit{10};

figure;
plotCurve({d1orig, d2orig}, splineData); 

%%
sdQuick = splineData;
sdQuick.N = splineData.N / 2;
sdQuick.varData.noPts = 3*sdQuick.N;
sdQuick.varData.pts = [];
sdQuick = constructKnots(sdQuick);
sdQuick = setupQuadData(sdQuick);

optQuick = struct('optDiff', true, ...
                  'optTra', true, ...
                  'optRot', true, ...
                  'optScal', true, ...
                  'useVarifold', true, ...
                  'varLambda', 100, ...
                  'useMultigrid', false, ...
                  'useAugmentedLagrangian', false, ...
                  'hansoNvec', 500, ...
                  'hansoMaxIt', 50, ...
                  'hansoCpuMax', 6, ...
                  'hansoNormTol', 1e-2, ...
                  'hansoPrtLevel', 0 );
              
options = struct( 'optDiff', true, ...
                  'optTra', true, ...
                  'optRot', true, ...
                  'optScal', true, ...
                  'useVarifold', true, ...
                  'varLambda', [], ...
                  'useMultigrid', true, ...
                  'useAugmentedLagrangian', false, ...
                  'hansoNvec', 500, ...
                  'hansoMaxIt', 500, ...
                  'hansoCpuMax', 60, ...
                  'hansoNormTol', 5e-2, ...
                  'hansoPrtLevel', 2 );
mgOptions = options;
mgOptions.useMultigrid = false;
mgOptions.hansoNormTol = 1e-1;
mgOptions.hansoCpuMax = 30;
options.mgOptions = mgOptions;

d1List = {};
d2List = {};
optEList = {};
optPathList = {};
optGaList = {};
splineDataList = {};
optionsList = {};
lambdaList = {};

%% Experiment setup
kernelSizeOrig = charLength / splineData.N * relKernelSize;

kList = {1, 2, 4, 8};
numExp = length(kList);

%% Do the computations
for jj = numExp:-1:1
    d1 = kList{jj} * d1orig;
    d2 = kList{jj} * d2orig;
    d1List{jj} = d1;
    d2List{jj} = d2;
    
    splineData.varData.kernelSizeGeom = kList{jj} * kernelSizeOrig;
    sdQuick.varData.kernelSizeGeom = splineData.varData.kernelSizeGeom;
    
    d1q = curveSpline2Spline(d1, splineData, sdQuick);
    d2q = curveSpline2Spline(d2, splineData, sdQuick);
    dInit = linearPath(d1q, d2q, sdQuick);
    
    [distEnEst, ~, ~, ~] = geodesicBvp(d1q, d2q, sdQuick, ...
                                       optQuick, 'initPath', dInit);
    
    % My estimate for lambda, assuming both parts of energy balance at the
    % optimum.
    lambda = distEnEst / varEnGoal;
    disp([num2str(distEnEst), ' ', num2str(sqrt(lambda))]);
    
    options.varLambda = lambda;
    options.mgOptions.varLambda = lambda;
    dInit = linearPath(d1, d2, splineData);
    lambdaList{jj} = lambda;
    
    [optE, optPath, optGa, ~] = geodesicBvp(d1, d2, splineData, ...
                                          options, 'initPath', dInit);
    optEList{jj} = optE;
    optPathList{jj} = optPath;
    optGaList{jj} = optGa;
    splineDataList{jj} = splineData;
    optionsList{jj} = options;

    analyzeBvpResults(d1, d2, optPath, splineData, options);
                                      
    %videoPath(['bump2_', num2str(jj), '.avi'], optPath, ...
    %          splineData, 'lineStyle', 'colour', 'lineWidth', 2);
end

%% Summarize results again
disp('------------ The Results ------------');
for jj = 1:numExp
    d1 = d1List{jj};
    d2 = d2List{jj};
    splineData = splineDataList{jj};
    optPath = optPathList{jj};
    optGa = optGaList{jj};
    
    analyzeBvpResults(d1, d2, optPath, ...
        splineData, optionsList{jj});
    
    dEnd = evalPath(optPath, 1, splineData);
    d2Ga = curveApplyGamma(d2, optGa, splineData);
    
    figure
    plotCurve({d2Ga, dEnd}, splineData, 'lineStyle', {'-k', '--r'});
end



