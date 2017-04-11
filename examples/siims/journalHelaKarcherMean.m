disp(mfilename);

%% Common parameters
splineData = constructEmptySplineData;
splineData.N = 40;
splineData.nS = 3;
splineData.Nt = 20;
splineData.nT = 2;
splineData.Nphi = 20;
splineData.nPhi = 3;
splineData.quadDegree = [6, 4];
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);
splineData.stepsT = 20;

options = struct( 'optDiff', true, ...
                  'optTra', true, ...
                  'optRot', true, ...
                  'optShift', true, ...
                  'tolFun', 1e-12, ...
                  'tolX', 1e-12, ...
                  'display', 'off', ... % 'off', 'iter-detailed'
                  'maxIter', 300, ...
                  'karcherTolGradNorm', 1e-3, ...
                  'karcherMaxIter', 50 );

dList = loadDataSet('hela_murphy', splineData, dataDir);
%indList = [1, 2, 3, 4, 5];
indList = 1:87;
initInd = 5;

%% Metric 1 - param
workDir = [prefixDir, 'hela_mean_1p/'];
resultsFile = [workDir, 'results.mat'];

if ~exist(workDir, 'dir')
    mkdir(workDir);
end

options.optDiff = false;
splineData.a = [1 0 0.125^4];

% Main calculation
journalHelaCalcKarcherMean(initInd, indList, splineData, options, ...
                           dataDir, workDir, resultsFile);

%% Metric 2 - param
workDir = [prefixDir, 'hela_mean_2p/'];
resultsFile = [workDir, 'results.mat'];

if ~exist(workDir, 'dir')
    mkdir(workDir);
end

options.optDiff = false;
splineData.a = [1 0 0.25^4];

% Main calculation
journalHelaCalcKarcherMean(initInd, indList, splineData, options, ...
                           dataDir, workDir, resultsFile);
  
%% Metric 3 - param
workDir = [prefixDir, 'hela_mean_3p/'];
resultsFile = [workDir, 'results.mat'];

if ~exist(workDir, 'dir')
    mkdir(workDir);
end

options.optDiff = false;
splineData.a = [1 0 0.5^4];

% Main calculation
journalHelaCalcKarcherMean(initInd, indList, splineData, options, ...
                           dataDir, workDir, resultsFile);

  
%% Metric 4 - param
workDir = [prefixDir, 'hela_mean_4p/'];
resultsFile = [workDir, 'results.mat'];

if ~exist(workDir, 'dir')
    mkdir(workDir);
end

options.optDiff = false;
splineData.a = [1 0 1^4];

% Main calculation
journalHelaCalcKarcherMean(initInd, indList, splineData, options, ...
                           dataDir, workDir, resultsFile);

%% Metric 1 - mod Diff
workDir = [prefixDir, 'hela_mean_1d/'];
resultsFile = [workDir, 'results.mat'];

if ~exist(workDir, 'dir')
    mkdir(workDir);
end

options.optDiff = true;
splineData.a = [1 0 0.125^4];

% Main calculation
journalHelaCalcKarcherMean(initInd, indList, splineData, options, ...
                           dataDir, workDir, resultsFile);
                      
%% Metric 2 - mod Diff
workDir = [prefixDir, 'hela_mean_2d/'];
resultsFile = [workDir, 'results.mat'];

if ~exist(workDir, 'dir')
    mkdir(workDir);
end

options.optDiff = true;
splineData.a = [1 0 0.25^4];

% Main calculation
journalHelaCalcKarcherMean(initInd, indList, splineData, options, ...
                           dataDir, workDir, resultsFile);
  
%% Metric 3 - mod Diff
workDir = [prefixDir, 'hela_mean_3d/'];
resultsFile = [workDir, 'results.mat'];

if ~exist(workDir, 'dir')
    mkdir(workDir);
end

options.optDiff = true;
splineData.a = [1 0 0.5^4];

% Main calculation
journalHelaCalcKarcherMean(initInd, indList, splineData, options, ...
                           dataDir, workDir, resultsFile);
  
%% Metric 4 - mod Diff
workDir = [prefixDir, 'hela_mean_4d/'];
resultsFile = [workDir, 'results.mat'];

if ~exist(workDir, 'dir')
    mkdir(workDir);
end

options.optDiff = true;
splineData.a = [1 0 1^4];

% Main calculation
journalHelaCalcKarcherMean(initInd, indList, splineData, options, ...
                           dataDir, workDir, resultsFile);