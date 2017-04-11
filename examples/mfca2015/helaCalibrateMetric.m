%% Prerequisites
% splineData, quadData
% prefixDir, splineCspDir, calibrateDir, calibrateFile

disp('Start calibration of constants in the metric.');

%% Load curves
prefix = prefixDir;
loaddir = splineCspDir;
calibrationfile = calibrateFile;

list_all = dir([prefix, loaddir, 'cell_*.mat']);
noCurves = length(list_all);

if ~isempty(maxNoCurves)
    noCurves = maxNoCurves;
end

%% Directory setup
if ~exist([prefix, savedir], 'dir')
    mkdir([prefix, savedir]);
end

%% Load all curves
dList = {};
noCurvesAll = noCurves;

for kk = noCurves:-1:1
    [d, ~] = loadCurve(list_all(kk).name, 'workdir', [prefix, loaddir]);
    dList{kk} = d;
end

%% Find average Riemannian energies
options = struct( 'optTra', true, ...
                  'optRot', true, ...
                  'optShift', true, ...
                  'rigidUseComp', false, ...
                  'rigidDisplay', 'off' );
                  
disp('Calculating Riemannian energies');
rlinEnergy = zeros([noCurves, noCurves, 4]);
for jj = noCurves:-1:2
    disp(jj);
    for kk = jj-1:-1:1
        [cAligned, ~] = rigidAlignment({dList{jj}, dList{kk}}, ...
            splineData, quadData, 'options', options);
        dLinearPath = linearPath(cAligned{1}, cAligned{2}, splineData);
        
        [~, comp] = pathRiemH2Energy( dLinearPath, splineData, ...
            quadData, quadDataTensor, 'a', [1 1 1] );
        rlinEnergy(kk,jj,3) = comp(3);
        rlinEnergy(kk,jj,2) = comp(2);
        rlinEnergy(kk,jj,1) = comp(1);
    end
end

for jj = 3:-1:1
    rlinMean(jj) = mean(nonzeros(rlinEnergy(:,:,jj)));
end

%% Perform calibration of metric constants
% Average L2, H1, H2 terms in the energy
meanL2 = rlinMean(1);
meanH1 = rlinMean(2);
meanH2 = rlinMean(3);

% With this scale L2 and H2 terms will be of equal size
lengthScale = (meanH2 / meanL2)^0.25;

% After we change the length scale, mean energies will scale too
meanL2 = meanL2 * lengthScale^3;
meanH1 = meanH1 * lengthScale;
meanH2 = meanH2 / lengthScale;

% With these coefficients all terms will be of equal size
coeffL2 = 1;
coeffH1 = sqrt(meanL2*meanH2) / meanH1;
coeffH2 = 1;

% Apply scaling to typical energy size
meanE = coeffL2 * meanL2 + coeffH1 * meanH1 + coeffH2 * meanH2;
coeffL2 = coeffL2 * typicalEnergy / meanE;
coeffH1 = coeffH1 * typicalEnergy / meanE;
coeffH2 = coeffH2 * typicalEnergy / meanE;

% Set coefficients in metric
aEqual = [coeffL2, coeffH1, coeffH2];
splineData.a = aEqual;

%% Calculate average length of curves
averageLength = 0;
for jj = noCurves:-1:1
    averageLength = averageLength + ...
        curveLength(dList{jj}, splineData, quadData);
end
averageLength = averageLength / noCurves;
averageLength = averageLength * lengthScale; % Length after rescaling

%% Save the found parameters
filename = [prefix, calibrationfile];
save(filename, 'aEqual', 'lengthScale', 'typicalEnergy', 'rlinMean', ...
               'averageLength');
disp(['typicalEnergy: ', num2str(typicalEnergy)]);
disp(['averageLength: ', num2str(averageLength)]);
disp(['Coefficients for equal weights: ', num2str(coeffL2), ' ', ...
    num2str(coeffH1), ' ', num2str(coeffH2)]);
disp('Calibration finished.');