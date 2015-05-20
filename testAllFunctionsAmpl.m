% Script for test testing all functionality
%% constructEmptySplineData, and setup parameters
splineData = constructEmptySplineData;
splineData.N = 20; %no. control points, must be bigger than n+1
splineData.Nt = 10 + 2; %Number of time control points
splineData.Nphi = 6; %No. control points for diffeomorphisms
splineData.nS = 4; %spacial degree
splineData.nT = 2; %time degree
splineData.nPhi = 3; %diffemorphism degree
splineData.quadDegree = [8,4]; %Quadrature precission
splineData.dSpace = 2;
splineData.a = [1 0 1];

%% constructKnots
splineData = constructKnots(splineData);

%% setupQuadData
tic
[quadData, quadDataTensor] = setupQuadData(splineData);
toc
%% constructSplineApproximation
f0 = @(t) [2*cos(t),sin(t)];
%f1 = @(t) 1/100*[150*cos(-t - 3*pi/4), 150*sin(-t - 3*pi/4)] ;
f1 = @(t) [cos(t),2*sin(t)+3];
f2 = @(t) [cos(t)+3,sin(t)-2];

d0 = constructSplineApproximation(f0,splineData);
d1 = constructSplineApproximation(f1,splineData);
d2 = constructSplineApproximation(f2,splineData);
x= pi/4;
A = [cos(x) sin(x);-sin(x) cos(x)];
d3 = d0;
for i=1:20
     d3(i,:) = (A*d0(i,:)')';
end


%% geodesicBVP
%[E_geo, dPath_optimal] = geodesicBvpAmpl(d0,d3,splineData,quadData,quadDataTensor,'datfileexists',true,'minrot',true);

%% Plot Path
%plotPath(dPath_optimal,splineData,8);

%% linearPath
%d_linear = linearPath(d0,d1,splineData);


%% geodesicForward (TODO)
%Nsteps = 100;
%v0 = pathVelocity(d_linear,0, splineData);
%q = geodesicForward(d0,d0,Nsteps,splineData,quadData);


%% Plotcurve
%d_array = {d0,d1};
% %[d_karcher,vel,grad_norm] = karcherMeanAmpl(d_array,splineData,quadData,quadDataTensor);
% plotCurve(d0,splineData);
% x = pi/2;
% A = [cos(x) -sin(x);sin(x) cos(x)];
% d3=d0;,
% for i=1:20
%     d3(i,:) = (A*d0(i,:)')'  
% end
% plotCurve(d3,splineData);


d_array = {d0,d1};
Dist = computeDistanceMatrix(d_array,splineData,quadData,quadDataTensor)



