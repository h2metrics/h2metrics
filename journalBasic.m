disp(mfilename);

splineDataList = {};
optionsList = {};
d1List = {};
d2List = {};
descriptionList = {};

% Uncomment for stand alone running
% dataDir = '~/diss/openprojects/h2_numerics/data/';
% prefixDir = '~/diss/openprojects/h2_numerics/journal/';
resultsFile = [prefixDir, 'basic_shapes.mat'];

%% Create directories
if ~exist(prefixDir, 'dir')
    mkdir(prefixDir);
end

%% Circle to wrap
splineData = constructEmptySplineData;
splineData.N = 40;
splineData.nS = 4;
splineData.Nt = 15;
splineData.nT = 3;
splineData.quadDegree = [6, 4];
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);
splineData.a = [1 0 1e-4];

d1 = loadDataSet('basic', splineData, '', ...
                 'curves', 'circle');
d2 = loadDataSet('basic', splineData, '', ...
                 'curves', 'wrap');
             
options = struct( 'optDiff', false, ...
                  'optTra', false, ...
                  'optRot', false, ...
                  'optShift', false, ...
                  'tolFun', 1e-12, ...
                  'tolX', 1e-12, ...
                  'display', 'off', ... % 'off', 'iter-detailed'
                  'maxIter', 300 );
              
splineDataList{end+1} = splineData;
optionsList{end+1} = options;
d1List{end+1} = d1;
d2List{end+1} = d2;
descriptionList{end+1} = 'Circle to wrap';

%% Propeller shapes
splineData = constructEmptySplineData;
splineData.N = 60;
splineData.nS = 4;
splineData.Nt = 15;
splineData.nT = 3;
splineData.quadDegree = [6, 4];
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);
splineData.a = [1 0 1e-4];

d1 = loadDataSet('basic', splineData, '', ...
                 'curves', 'prop3', 'noise', 0.);
d2 = loadDataSet('basic', splineData, '', ...
                 'curves', 'prop4', 'noise', 0.);

options = struct( 'optDiff', false, ...
                  'optTra', false, ...
                  'optRot', false, ...
                  'optShift', false, ...
                  'tolFun', 1e-12, ...
                  'tolX', 1e-12, ...
                  'display', 'off', ... % 'off', 'iter-detailed'
                  'maxIter', 300 );
              
splineDataList{end+1} = splineData;
optionsList{end+1} = options;
d1List{end+1} = d1;
d2List{end+1} = d2;
descriptionList{end+1} = 'Propeller 3 to 4';

%% Corpus Callosum
splineData = constructEmptySplineData;
splineData.N = 60; %N=100
splineData.nS = 4;
splineData.Nt = 15;
splineData.nT = 3;
splineData.quadDegree = [6, 4];
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);
splineData.a = [1 0 1e-5];

rigidOptions = struct( 'optDiff', false, ...
                       'optTra', false, ...
                       'optRot', false, ...
                       'optShift', true );

dList = loadDataSet( 'corpus_callosum_tilak', splineData, dataDir, ...
                     'class', 'normal', 'ind', [1, 2], 'constSpeed');
d1 = dList{1}; d2 = dList{2};
meanLength = ( curveLength(d1, splineData) + ...
               curveLength(d2, splineData) ) / 2;
d1 = d1 / meanLength * 2*pi;
d2 = d2 / meanLength * 2*pi;
[dList, ~] = rigidAlignment({d1, d2}, splineData, 'options', rigidOptions);
d2 = dList{2};

options = struct( 'optDiff', false, ...
                  'optTra', true, ...
                  'optRot', true, ...
                  'optShift', true, ...
                  'tolFun', 1e-12, ...
                  'tolX', 1e-12, ...
                  'display', 'off', ... % 'off', 'iter-detailed'
                  'maxIter', 300, ...
                  'rigidGlobalRot', false, ...
                  'rigidUseComp', true );
              
splineDataList{end+1} = splineData;
optionsList{end+1} = options;
d1List{end+1} = d1;
d2List{end+1} = d2;
descriptionList{end+1} = 'Corpus Callosum';

%% HeLa cells
splineData = constructEmptySplineData;
splineData.N = 40; %N=100
splineData.nS = 3;
splineData.Nt = 20;
splineData.nT = 2;
splineData.Nphi = 20;
splineData.nPhi = 3;
splineData.quadDegree = [6, 4];
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);
splineData.a = [1 0 2^-12];
splineData.stepsT = 20;

rigidOptions = struct( 'optDiff', false, ...
                       'optTra', false, ...
                       'optRot', false, ...
                       'optShift', true );

dList = loadDataSet( 'hela_murphy', splineData, dataDir, 'ind', [3, 24]);
d1 = dList{1}; d2 = dList{2};
meanLength = ( curveLength(d1, splineData) + ...
               curveLength(d2, splineData) ) / 2;
d1 = d1 / meanLength * 2*pi;
d2 = d2 / meanLength * 2*pi;
[dList, ~] = rigidAlignment({d1, d2}, splineData, 'options', rigidOptions);
d2 = dList{2};

options = struct( 'optDiff', true, ...
                  'optTra', true, ...
                  'optRot', true, ...
                  'optShift', true, ...
                  'tolFun', 1e-12, ...
                  'tolX', 1e-12, ...
                  'display', 'off', ... % 'off', 'iter-detailed'
                  'maxIter', 300, ...
                  'rigidGlobalRot', false, ...
                  'rigidUseComp', true );
              
splineDataList{end+1} = splineData;
optionsList{end+1} = options;
d1List{end+1} = d1;
d2List{end+1} = d2;
descriptionList{end+1} = 'HeLa cells';

%% Save the data
% Strip quadData
for jj = 1:length(splineDataList)
    splineDataList{jj}.quadData = [];
    splineDataList{jj}.quadDataTensor = [];
end
save(resultsFile, 'splineDataList', 'd1List', 'd2List', ...
        'descriptionList', 'optionsList');

%% Make some basic computations
noExamples = length(splineDataList);
for jj = 1:noExamples
    splineData = splineDataList{jj};
    splineData = setupQuadData(splineData);
    quadData = splineData.quadData;
    quadDataTensor = splineData.quadDataTensor;
    options = optionsList{jj};
    d1 = d1List{jj};
    d2 = d2List{jj};
    
    % Describe splineData
    disp(' ');
    disp('-----------------------------------------------');
    disp(['Example ', num2str(jj), ' - ', descriptionList{jj}]);
    disp(['N=', num2str(splineData.N), ', nS=', num2str(splineData.nS)]);
    disp(['Nt=', num2str(splineData.Nt), ', nT=', num2str(splineData.nT)]);
    disp(['Nphi=', num2str(splineData.Nphi), ...
          ', nPhi=', num2str(splineData.nPhi)]);
    quadPoints = ceil((splineData.quadDegree-1)/2);
    disp(['quadDegree=(', num2str(splineData.quadDegree(1)), ', ',...
          num2str(splineData.quadDegree(2)), '), ',...
          'quadPoints=(', num2str(quadPoints(1)), ', ', ...
          num2str(quadPoints(2)), ')']);
    disp(['a(1)=', num2str(splineData.a(1)), ...
          ', a(2)=', num2str(splineData.a(2)), ...
          ', a(2)=', num2str(splineData.a(3))]);
    disp(['optDiff=', num2str(options.optDiff), ...
          ', optTra=', num2str(options.optTra), ...
          ', optRot=', num2str(options.optRot), ...
          ', optShift=', num2str(options.optShift) ]);
      
    % Curve information
    disp(' ');
    disp(['Curve 1, Len=', ...
          num2str(curveLength(d1, splineData))]);
    disp(['Curve 2, Len=', ...
          num2str(curveLength(d2, splineData))]);
      
    % Calculate geodesic
    tic;
    [optE, optPath, optGa, info] = geodesicBvp(d1, d2, splineData, ...
                                               'options', options);
    time = toc;
    dist = sqrt(optE);
    [E, comp] = pathRiemH2Energy(optPath, splineData);
    E_L2 = comp(1) / E;
    E_H1 = comp(2) / E;
    E_H2 = comp(3) / E;
    
    disp(' ');
    disp('Geodesic BVP');
    disp(['dist(d1, d2)=', num2str(dist), ', Energy=', num2str(optE)]);
    disp(['E_L2=', num2str(E_L2), ', E_H1=', num2str(E_H1), ...
          ', E_H2=', num2str(E_H2)]);
    disp(['Computation time: ', num2str(time)]);
    disp(['Iterations: ', num2str(info.noIter)]);
    
end