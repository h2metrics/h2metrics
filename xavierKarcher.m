% Load curves
load('xavier_splines_v1_N20','dPatients');
splineData = constructEmptySplinedata;
splineData.N = 20; %no. control points, must be bigger than n+1
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

[dPatients_karcher,vel,grad_norm] = karcherMeanAmpl(dPatients,splineData,quadData,quadDataTensor);

%Save results
save('xavierKarcher.mat','dPatients_karcher','vel','grad_norm');