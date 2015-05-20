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
[C,vel,grad_norm] = karcherMeanAmpl(dSick,splineData,quadData,quadDataTensor);
plotKarcherWithCurves(C,dHealthy,splineData)


