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
%A1B1C1
splineData.a = [A B C];
[CSick,vel,grad_norm] = karcherMeanAmpl(dSick,splineData,quadData,quadDataTensor);
save('CSickA1B1C1','CSick');
[CHealthy,vel,grad_norm] = karcherMeanAmpl(dHealthy,splineData,quadData,quadDataTensor);
save('CHealthyA1B1C1','CHealthy');

%A80B10C10
splineData.a = [80*A 10*B 10*C];
[CSick,vel,grad_norm] = karcherMeanAmpl(dSick,splineData,quadData,quadDataTensor);
save('CSickA80B10C10','CSick');
[CHealthy,vel,grad_norm] = karcherMeanAmpl(dHealthy,splineData,quadData,quadDataTensor);
save('CHealthyA80B10C10','CHealthy');


%% Calculate Geodesic between Karchermeans
%A1B1C1
splineData.a = [A B C];
load('CHealthyA1B1C1');
load('CSickA1B1C1');
[~,dPath] = geodesicBvpAmpl(CSick,CHealthy,splineData,quadData,quadDataTensor);
save('dPathA1B1C1MeanToMean','dPath');

%A80B10C10
splineData.a = [80*A 10*B 10*C];
load('CHealthyA80B10C10');
load('CSickA80B10C10');
[~,dPath] = geodesicBvpAmpl(CSick,CHealthy,splineData,quadData,quadDataTensor);
save('dPathA80B10C10MeanToMean','dPath');


%% Compute distance Matrix with Karcher means
%A1B1C1
splineData.a = [A B C];
load('CHealthyA1B1C1');
load('CSickA1B1C1');
dPatientsMeans = {CHealthy,dPatients_noloops{1:19},CSick};
DistA1B1C1Means = computeDistanceMatrix(dPatientsMeans,splineData,quadData,quadDataTensor);
save('DistA1B1C1Means','DistA1B1C1');

%A80B10C10
splineData.a = [80*A 10*B 10*C];
load('CHealthyA80B10C10');
load('CSickA80B10C10');
dPatientsMeans = {CHealthy,dPatients_noloops{1:19},CSick};
DistA1B1C1Means = computeDistanceMatrix(dPatientsMeans,splineData,quadData,quadDataTensor);
save('DistA1B1C1Means','DistA1B1C1');


















