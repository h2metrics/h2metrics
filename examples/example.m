clear;

%% set up spline interpolation
splineData = constructEmptySplineData;
splineData.N = 40;          % Spatial control points
splineData.nS = 4;          % Spatial spline degree
splineData.Nt = 15;         % Temporal control points
splineData.nT = 3;          % Temporal spline degree
splineData.Nphi = 20;       % Control points for diffeomorphism
splineData.nPhi = 3;        % Spline degree for diffeomorphism
splineData.quadDegree = [6, 4]; % Quadratue orders in space and time
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);
splineData.a = [1 0 1e-3];  % Constants in the metric
splineData.stepsT = 10;     % Time steps for initial value problem

%% Load a selection of sample curves
dList = loadDataSetBasic(splineData, 'curves', {'prop3', 'prop4'});
d1 = dList{1};  % Initial curve
d2 = dList{2};  % Final curve

%% Plot the boundary curves
figure();
plotCurve(d1, splineData);
plotCurve(d2, splineData);

%% Solve a geodesic boundary value problem between parametrized curves
% Note: neither reparametrizations nor Euclidean motions are factored out
[optE, ~, ~, info] = geodesicBvp(d1, d2, splineData);

disp(['Number of iterations: ', num2str(info.noIter)]);
disp(['Geodesic distance: ', num2str(sqrt(optE))]);

%% Customize options using optional arguments to geodesicBvp
options = struct( ...
    'optDiff', true, ...        % Reparameterizations
    'optTra', true, ...         % Translations
    'optRot', true, ...         % Rotations
    'optShift', true, ...       % Shifts of the parametrization
    'tolFun', 1e-12, ...        % Tolerance for optimization routine
    'tolX', 1e-12, ...          % Tolerance for optimization routine
    'display', 'iter-detailed', ... % Use 'off' or 'iter-detailed'
    'maxIter', 300 );           % Max. # of iterations
              
%% Solve geodesic BVP again
[optE, optPath, optGa, info] = geodesicBvp(d1, d2, splineData, ...
                                           'options', options);

%% Plot the geodesic
figure();
plotPath(optPath, splineData);

%% Compute Karcher mean
% This requires the Manopt library from http://www.manopt.org/

options.optDiff = false;
options.optTra = false;
options.optRot  = false; 
options.optShift = false;
options.display = 'off';

options.karcherTolGradNorm = 1e-3;  % Tolerance for conjugate gradients
options.karcherMaxIter = 10;        % Max. # of iterations for CG

dList = loadDataSetBasic(splineData, 'curves', {'prop3', 'prop4'});
if exist('conjugategradient','file') == 2 % Only exectute if Manopt exists
    [dMean, info] = karcherMeanManopt( dList, ...
        splineData, 'options', options, 'meanInit', dList{1} );
end

%% Plot the curves and the Karcher mean
if exist('dMean', 'var') == 1
    figure();
    plotCurve(dList{1}, splineData, 'k--');
    plotCurve(dList{2}, splineData, 'k--');
    plotCurve(dMean, splineData, 'b-');
end