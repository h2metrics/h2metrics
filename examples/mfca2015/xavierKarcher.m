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
%%
% splineData.a = [10*A 90*B 0];
% [CSick,vel,grad_norm] = karcherMeanAmpl(dSick,splineData,quadData,quadDataTensor);
% save('CSickA10B90C00','CSick');
plotKarcherWithCurves(CSick,dSick,splineData)
export_fig CSickA10B90C00.pdf

[CHealthy,vel,grad_norm] = karcherMeanAmpl(dHealthy,splineData,quadData,quadDataTensor);
save('CHealthyA10B90C00','CHealthy');
plotKarcherWithCurves(CHealthy,dHealthy,splineData)
export_fig CHealthyA10B90C00.pdf
% 
% %%
% splineData.a = [80*A 10*B 10*C];
% [CSick,vel,grad_norm] = karcherMeanAmpl(dSick,splineData,quadData,quadDataTensor);
% save('CSickA80B10C10','CSick');
% plotKarcherWithCurves(CSick,dSick,splineData)
% export_fig CSickA80B10C10.pdf
% 
% [CHealthy,vel,grad_norm] = karcherMeanAmpl(dHealthy,splineData,quadData,quadDataTensor);
% save('CHealthyA80B10C10','CHealthy');
% plotKarcherWithCurves(CHealthy,dHealthy,splineData)
% export_fig CHealthyA80B10C10.pdf
% 
% %%
% splineData.a = [A B 98*C];
% [CSick,vel,grad_norm] = karcherMeanAmpl(dSick,splineData,quadData,quadDataTensor);
% save('CSickA1B1C98','CSick');
% plotKarcherWithCurves(CSick,dSick,splineData)
% export_fig CSickA1B1C98.pdf
% 
% [CHealthy,vel,grad_norm] = karcherMeanAmpl(dHealthy,splineData,quadData,quadDataTensor);
% save('CHealthyA1B1C98','CHealthy');
% plotKarcherWithCurves(CHealthy,dHealthy,splineData)
% export_fig CHealthyA1B1C98.pdf

%%Calculate distances to the Karcher mean
% load('CHealthyA1B1C1');
% Dist = matchOneToAll(CHealthy,dHealthy,splineData,quadData,quadDataTensor);
% save('DistA1B1C1dHealthy','Dist');
% 
% load('CSickA1B1C1');
% Dist = matchOneToAll(CSick,dSick,splineData,quadData,quadDataTensor);
% save('DistA1B1C1dSick','Dist');
% load('DistA1B1C1dSick');
% mean(Dist)
% std(Dist)
% var(Dist)
% 
% 
% load('CHealthyA1B1C1');
% Dist = matchOneToAll(CHealthy,dHealthy,splineData,quadData,quadDataTensor);
% save('DistA1B1C1dHealthy','Dist');
% load('DistA1B1C1dHealthy');
% mean(Dist)
% std(Dist)
% var(Dist)


%%Plot Geodesic between Karchermeans
% load('CHealthyA1B1C1');
% load('CSickA1B1C1');
% [~,dPath] = geodesicBvpAmpl(CSick,CHealthy,splineData,quadData,quadDataTensor);
% 
% plotParticlePath(dPath,splineData,6,40);
% export_fig MeanSMeanHA1B1C1.pdf
% plotPathSingle(dPath,splineData,6)
% plotPathRow(dPath,splineData,6);
% export_fig MeanSMeanHA1B1C1Row.pdf










