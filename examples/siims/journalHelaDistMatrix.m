disp(mfilename);
clear;

dataDir = '~/diss/openprojects/h2_numerics/data/';
prefixDir = '~/diss/openprojects/h2_numerics/journal/';

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

% indList = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
indList = 1:87;

meanInitInd = 5;

%% Create curves
noCurves = length(indList);                
dList = loadDataSet('hela_murphy', splineData, dataDir, 'ind', indList);
                
% Center curves and normalize average length
avLen = 0;
for jj = noCurves:-1:1
    avLen = avLen + curveLength(dList{jj}, splineData);
end
avLen = avLen / noCurves;

for jj = noCurves:-1:1
    dList{jj} = curveCenter(dList{jj}, splineData);
    dList{jj} = dList{jj} / avLen * 2*pi;
end

%% Metric 1 - mod Diff
results = matfile([prefixDir, 'hela_new/hela_1d.mat'], 'Writable', true);

options.optDiff = true;
splineData.a = [1 0 2^-12];

results.splineData = splineData;
results.options = options;
results.dList = dList;
results.indList = indList;
results.noCurves = length(indList);

results.meanInitInd = meanInitInd;

journalHelaCalcDistMatrix(results, true);
journalHelaCalcKarcherMean(results);

%% Metric 2 - mod Diff
results = matfile([prefixDir, 'hela_new/hela_2d.mat'], 'Writable', true);

options.optDiff = true;
splineData.a = [1 0 2^-8];

results.splineData = splineData;
results.options = options;
results.dList = dList;
results.indList = indList;
results.noCurves = length(indList);

results.meanInitInd = meanInitInd;

journalHelaCalcDistMatrix(results, true);
journalHelaCalcKarcherMean(results);


% %% Metric 1 - param
% workDir = [prefixDir, 'hela_dist_1p/'];
% resultsFile = [workDir, 'results.mat'];
% 
% if ~exist(workDir, 'dir')
%     mkdir(workDir);
% end
% 
% options.optDiff = false;
% splineData.a = [1 0 0.125^4];
% 
% % Main calculation
% journalHelaCalcDistMatrix(indList, splineData, options, ...
%                           dataDir, workDir, resultsFile, false);
% 
% %% Metric 2 - param
% workDir = [prefixDir, 'hela_dist_2p/'];
% resultsFile = [workDir, 'results.mat'];
% 
% if ~exist(workDir, 'dir')
%     mkdir(workDir);
% end
% 
% options.optDiff = false;
% splineData.a = [1 0 0.25^4];
% 
% % Main calculation
% journalHelaCalcDistMatrix(indList, splineData, options, ...
%                           dataDir, workDir, resultsFile, false);
%   
% % %% Metric 3 - param
% % workDir = [prefixDir, 'hela_dist_3p/'];
% % resultsFile = [workDir, 'results.mat'];
% % 
% % if ~exist(workDir, 'dir')
% %     mkdir(workDir);
% % end
% % 
% % options.optDiff = false;
% % splineData.a = [1 0 0.5^4];
% % 
% % % Main calculation
% % journalHelaCalcDistMatrix(indList, splineData, options, ...
% %                           dataDir, workDir, resultsFile, false);
% %   
% % %% Metric 4 - param
% % workDir = [prefixDir, 'hela_dist_4p/'];
% % resultsFile = [workDir, 'results.mat'];
% % 
% % if ~exist(workDir, 'dir')
% %     mkdir(workDir);
% % end
% % 
% % options.optDiff = false;
% % splineData.a = [1 0 1^4];
% % 
% % % Main calculation
% % journalHelaCalcDistMatrix(indList, splineData, options, ...
% %                           dataDir, workDir, resultsFile, false);
% 
% %% Metric 1 - mod Diff
% workDir = [prefixDir, 'hela_dist_1d/'];
% resultsFile = [workDir, 'results.mat'];
% 
% if ~exist(workDir, 'dir')
%     mkdir(workDir);
% end
% 
% options.optDiff = true;
% splineData.a = [1 0 0.125^4];
% 
% % Main calculation
% journalHelaCalcDistMatrix(indList, splineData, options, ...
%                           dataDir, workDir, resultsFile, true);
%                       
% %% Metric 2 - mod Diff
% workDir = [prefixDir, 'hela_dist_2d/'];
% resultsFile = [workDir, 'results.mat'];
% 
% if ~exist(workDir, 'dir')
%     mkdir(workDir);
% end
% 
% options.optDiff = true;
% splineData.a = [1 0 0.25^4];
% 
% % Main calculation
% journalHelaCalcDistMatrix(indList, splineData, options, ...
%                           dataDir, workDir, resultsFile, true);
%   
% % %% Metric 3 - mod Diff
% % workDir = [prefixDir, 'hela_dist_3d/'];
% % resultsFile = [workDir, 'results.mat'];
% % 
% % if ~exist(workDir, 'dir')
% %     mkdir(workDir);
% % end
% % 
% % options.optDiff = true;
% % splineData.a = [1 0 0.5^4];
% % 
% % % Main calculation
% % journalHelaCalcDistMatrix(indList, splineData, options, ...
% %                           dataDir, workDir, resultsFile, true);
% %   
% % %% Metric 4 - mod Diff
% % workDir = [prefixDir, 'hela_dist_4d/'];
% % resultsFile = [workDir, 'results.mat'];
% % 
% % if ~exist(workDir, 'dir')
% %     mkdir(workDir);
% % end
% % 
% % options.optDiff = true;
% % splineData.a = [1 0 1^4];
% % 
% % % Main calculation
% % journalHelaCalcDistMatrix(indList, splineData, options, ...
% %                           dataDir, workDir, resultsFile, true);