%% Prerequisites
% splineData, quadData
% weightL2, weightH1, weightH2
% prefixDir, calibrateFile, splineCspDir, paramDir

%% Load curves
prefix = prefixDir;
calibrationfile = calibrateFile;
loaddir = splineCspDir;
savedir = [paramDir, 'spline/'];

list_all = dir([prefix, loaddir, 'cell_*.mat']);
noCurves = length(list_all);

if ~isempty(maxNoCurves)
    noCurves = maxNoCurves;
end
                      
%% Load calibration file
calibrationMatFile = matfile([prefix, calibrationfile]);
aEqual = calibrationMatFile.aEqual;
lengthScale = calibrationMatFile.lengthScale;

% Apply weights
coeffL2 = weightL2 * aEqual(1);
coeffH1 = weightH1 * aEqual(2);
coeffH2 = weightH2 * aEqual(3);

% Set coefficients in metric
a = [coeffL2, coeffH1, coeffH2];
splineData.a = a;

%% Directory setup
if ~exist([prefix, savedir], 'dir')
    mkdir([prefix, savedir]);
end

%% Load all curves
dList = {};

for kk = noCurves:-1:1
    [d, ~] = loadCurve(list_all(kk).name, 'workdir', [prefix, loaddir]);
    dList{kk} = d;
end

%% Center curves
dAligned = {};
for kk = noCurves:-1:1
    dAligned{kk} = curveCenter(dList{kk}, splineData, quadData);
end

%% Apply length calibration
for kk = noCurves:-1:1
    dAligned{kk} = lengthScale * dAligned{kk};
end

%% Save curves again
for kk = 1:noCurves
    saveCurve(list_all(kk).name, dAligned{kk}, splineData, ...
        'workdir', [prefix, savedir], 'manualName', true);
end

disp('Calibration applied.');
