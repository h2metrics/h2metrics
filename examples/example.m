%% This source contains a verty basic example for the functionality of the code
%% Precaution
clear all;
%% Set paths
addpath(genpath('../source'));
addpath('../lib/varifolds');
addpath('../lib/hanso2_2');
addpath('../lib/export_fig');
%% Setup parameters for curve discretization
splineData = constructSplineData;
splineData.N = 40; % Number of spline controll points for spatial discretization
splineData.nS = 2; % Quadratic splines are needed for second order metrics  
splineData.Nt = 5; % Number of spline controll points for time discretization 
splineData.curveClosed=1; %Closed Curves (set to 0 for open curves)
splineData = finishSetup(splineData);
splineData.varData.noPts = 100; % Number of points for similarity measure (varifold) evaulation
splineData.varData.kernelSizeGeom = 0.1;
splineData.varData.kernelSizeGrass = 0.3;
splineData.varData.pts = [];
splineData = finishSetup( splineData );
%% Load some curves
c0 = importdata('./data/OAS1_0016.txt');
c1 = importdata('./data/OAS1_0022.txt');
c2 = importdata('./data/OAS1_0023.txt');
%Construct spline approximation
d0= constructSplineApproximation(c0(2:end,:),splineData);
d1= constructSplineApproximation(c1(2:end,:),splineData);
d2= constructSplineApproximation(c2(2:end,:),splineData);
%Rescale curves to length 2pi (not necessary)
d0 = 2*pi* d0/curveLength(d0,splineData);
d1 = 2*pi*d1/curveLength(d1,splineData);
d2= 2*pi*d2/curveLength(d2,splineData);
dList= {d0,d1,d2};
%Prealign the data (not necessary)
dList = rigidAlignmentVarifold({d0,d1,d2},splineData); 
plotCurve(dList, splineData) %Plot the data
%% Calculate Geodesic between c0 and c1
splineData.a=[0 1 0 0 0]; %Constants for the metric:
splineData.scaleInv=1; %length weighted metric (set to zero for constant coeff. metrics)
splineData.options.optDiff=1;  %Minimize over reparamtrizations of the target curve
splineData.options.optScal=1; %Minimize over scalings of the target curve
splineData.options.optRot=1; %Minimize over rotations of the target curve
splineData.options.optTra=1; %Minimize over translations of the target curve
splineData.options.useAugmentedLagrangian = false; %Use quadratic penalty term 
splineData.options.varLambda = 100; %Weight of the similarity measure
%% Minimize and plot optimal geodesic
tic
[optE, optPath, optGa, ~] = geodesicBvp(dList{1}, 1.1*dList{2}, ...
    splineData ,splineData.options);
toc;
dEnd=curveApplyGamma(dList{2}, optGa, splineData);
disp('Geodesic distance between d0 and d1 is:');
disp((optE)^(1/2));
% Plot of the minimizing geodesic. The target curve is ploted in red.
clf
plotPath2(optPath,dEnd,splineData)
%% Compute Karcher Mean (needs to setup manopt, 
%run first the importmanopt script in the lib/manopt folder
karcherOptions.karcherMaxIter = 10;
dMean = karcherMeanManopt(dList, splineData,'options',karcherOptions);
%% Plot mean in red and curves
plotCurve(dMean, splineData,'lineStyle', 'r-')
plotCurve(dList, splineData)
%% Tangent space PCA
[U,Lambda,G, vList] = TangentPCA(dList,dMean,splineData);
VisualizePCA(U,Lambda,G, vList,splineData)
