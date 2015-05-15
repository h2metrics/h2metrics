% Load curves
load('xavier_splines_v1_N30','dPatients');
splineData = constructEmptySplinedata;
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
dSick = {dPatients{16:25}};
splineData.a = [1 1 0];
% DistA1B1C0 = computeDistanceMatrix(dPatients,splineData,quadData,quadDataTensor);
% save('DistA1B1C0','DistA1B1C0');
% splineData.a = [1 0.1 0];
% DistA1B01C0 = computeDistanceMatrix(dPatients,splineData,quadData,quadDataTensor);
% save('DistA1B01C0','DistA1B01C0');
% splineData.a = [1 0.003 0];
% DistA1B003C0 = computeDistanceMatrix(dPatients,splineData,quadData,quadDataTensor);
% save('DistA1B003C0','DistA1B003C0');
% 






%[dPatients_karcher,vel,grad_norm] = karcherMeanAmpl(dHealthy,splineData,quadData,quadDataTensor);
plotCurves({dPatients{8}},splineData);
%[E_geo, dPath] = geodesicBvpAmpl(dPatients{1},dPatients{25},splineData,quadData,quadDataTensor,'datfileexists',true);
%plotPath(dPath,splineData,8);
 
%dPath = linearPath(dPatients{16},dPatients{20},splineData);
%[G, comp] = pathRiemH2Energy( dPath, splineData, quadData, quadDataTensor);
%[A B C] = determineConstants(dSick,splineData,quadData,quadDataTensor);