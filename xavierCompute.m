load('xavier_splines_v1_N30_noloops','dPatients_noloops');
splineData = constructEmptySplineData;
splineData.N = 30; %no. control points, must be bigger than n+1
splineData.Nt = 10 + 2; %Number of time control points
splineData.Nphi = 6; %No. control points for diffeomorphisms
splineData.nS = 3; %spacial degree
splineData.nT = 2; %time degree
splineData.nPhi = 3; %diffemorphism degree
splineData.quadDegree = [8,4]; %Quadrature precission
splineData.dSpace = 2;
splineData.noInterpolS = 2 * splineData.N; % For composition
splineData = constructKnots(splineData);
[quadData, quadDataTensor] = setupQuadData(splineData);
dSick = {dPatients_noloops{10:19}};
dHealthy = {dPatients_noloops{1:9}};
splineData.a = [1 1 1];
[A B C] = determineConstants(dPatients_noloops,splineData,quadData,quadDataTensor);
splineData.a = [A B C];


%% Calculate Karcher Means
%%A1B1C1
splineData.a = [A B C];
[CSick,vel,grad_norm] = karcherMeanAmpl(dSick,splineData,quadData,quadDataTensor);
save('CSickA1B1C1','CSick');
[CHealthy,vel,grad_norm] = karcherMeanAmpl(dHealthy,splineData,quadData,quadDataTensor);
save('CHealthyA1B1C1','CHealthy');
load('CHealthyA1B1C1');
load('CSickA1B1C1');
[CInitialGuess,vel,grad_norm] = karcherMeanAmpl({CSick,CHealthy},splineData,quadData,quadDataTensor);
save('CInitialGuess','CInitialGuess');
[CPatients,vel,grad_norm] = karcherMeanAmplGuess(dPatients_noloops,CInitialGuess,splineData,quadData,quadDataTensor);
save('CPatientsA1B1C1','CPatients');

%% Calculate Initial velocities from the mean to subjects
load('CHealthyA1B1C1')
[DistToMeanHealthyA1B1C1,InitialVelHealthyA1B1C1] = matchOneToAll(CHealthy,dHealthy,splineData,quadData,quadDataTensor);
save('DistToMeanHealthyA1B1C1','DistToMeanHealthyA1B1C1');

load('CSickA1B1C1')
[DistToMeanSickA1B1C1,InitialVelSickA1B1C1] = matchOneToAll(CSick,dSick,splineData,quadData,quadDataTensor);
save('DistToMeanSickA1B1C1','DistToMeanSickA1B1C1');

load('CPatientsA1B1C1')
[DistToMeanPatientsA1B1C1,InitialVelPatientsA1B1C1] = matchOneToAll(CPatients,dPatients_noloops,splineData,quadData,quadDataTensor);
save('DistToMeanPatientsA1B1C1','DistToMeanPatientsA1B1C1');
save('InitialVelPatientsA1B1C1','InitialVelPatientsA1B1C1');

%% Calculate Distance between means
Dist = matchOneToAll(CSick,{CHealthy},splineData,quadData,quadDataTensor);

%% Compute Distance-Matrix with Karcher means
splineData.a = [A B C];
load('CHealthyA1B1C1');
load('CSickA1B1C1');
dPatientsMeans = {CHealthy,dPatients_noloops{1:19},CSick};
DistA1B1C1Means = computeDistanceMatrix(dPatientsMeans,splineData,quadData,quadDataTensor);
save('DistA1B1C1Means','DistA1B1C1Means');


%% Other Values for the constants in the  metric
% %A50B30C20
% splineData.a = [50*A 30*B 20*C];
% [CSick,vel,grad_norm] = karcherMeanAmpl(dSick,splineData,quadData,quadDataTensor);
% save('CSickA50B30C20','CSick');
% [CHealthy,vel,grad_norm] = karcherMeanAmpl(dHealthy,splineData,quadData,quadDataTensor);
% save('CHealthyA50B30C20','CHealthy');
% %A10B10C80
% splineData.a = [10*A 10*B 80*C];
% [CSick,vel,grad_norm] = karcherMeanAmpl(dSick,splineData,quadData,quadDataTensor);
% save('CSickA10B10C80','CSick');
% [CHealthy,vel,grad_norm] = karcherMeanAmpl(dHealthy,splineData,quadData,quadDataTensor);
% save('CHealthyA10B10C80','CHealthy');
% %A20B30C50
% splineData.a = [20*A 30*B 50*C];
% [CSick,vel,grad_norm] = karcherMeanAmpl(dSick,splineData,quadData,quadDataTensor);
% save('CSickA20B30C50','CSick');
% [CHealthy,vel,grad_norm] = karcherMeanAmpl(dHealthy,splineData,quadData,quadDataTensor);
% save('CHealthyA20B30C50','CHealthy');
% 
% 
% %% Calculate Geodesic between Karchermeans
% %A1B1C1
% splineData.a = [A B C];
% load('CHealthyA1B1C1');
% load('CSickA1B1C1');
% [~,dPath] = geodesicBvpAmpl(CSick,CHealthy,splineData,quadData,quadDataTensor);   
% save('dPathA1B1C1MeanToMean','dPath');
% 
% %% Compute distance Matrix with Karcher means
% %A1B1C1
% %A50B30C20
% splineData.a = [50*A 30*B 20*C];
% load('CHealthyA50B30C20');
% load('CSickA50B30C20');
% dPatientsMeans = {CHealthy,dPatients_noloops{1:19},CSick};
% DistA50B30C20Means = computeDistanceMatrix(dPatientsMeans,splineData,quadData,quadDataTensor);
% save('DistA50B30C20Means','DistA50B30C20Means');
% %A10B10C80
% splineData.a = [10*A 10*B 80*C];
% load('CHealthyA10B10C80');
% load('CSickA10B10C80');
% dPatientsMeans = {CHealthy,dPatients_noloops{1:19},CSick};
% DistA10B10C80Means = computeDistanceMatrix(dPatientsMeans,splineData,quadData,quadDataTensor);
% save('DistA10B10C80Means','DistA10B10C80Means');
% %A20B30C50
% splineData.a = [20*A 30*B 50*C];
% load('CHealthyA20B30C50');
% load('CSickA20B30C50');
% dPatientsMeans = {CHealthy,dPatients_noloops{1:19},CSick};
% DistA20B30C50Means = computeDistanceMatrix(dPatientsMeans,splineData,quadData,quadDataTensor);
% save('DistA20B30C50Means','DistA20B30C50Means');
% 
















