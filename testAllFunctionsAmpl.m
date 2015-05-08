% Script for test testing all functionality
%% constructEmptySplineData, and setup parameters
splineData = constructEmptySplinedata;
splineData.N = 20; %no. control points, must be bigger than n+1
splineData.Nt = 10 + 2; %Number of time control points
splineData.Nphi = 6; %No. control points for diffeomorphisms
splineData.nS = 4; %spacial degree
splineData.nT = 2; %time degree
splineData.nPhi = 3; %diffemorphism degree
splineData.quadDegree = [8,4]; %Quadrature precission
splineData.dSpace = 2;

%% constructKnots
splineData = constructKnots(splineData);

%% setupQuadData
tic
[quadData, quadDataTensor] = setupQuadData(splineData);
toc
%% constructSplineApproximation
f0 = @(t) [cos(t),sin(t)];
%f1 = @(t) 1/100*[150*cos(-t - 3*pi/4), 150*sin(-t - 3*pi/4)] ;
f1 = @(t) [cos(t)+5,sin(t)+5];
f2 = @(t) [cos(t)+3,sin(t)-2];

d0 = constructSplineApproximation(f0,splineData);
d1 = constructSplineApproximation(f1,splineData);
d2 = constructSplineApproximation(f2,splineData);

%% geodesicBVP
[E_geo, dPath_optimal] = geodesicBvpAmpl(d0,d1,splineData,quadData,quadDataTensor,'datfileexists',false);

%% Plot Path
%plotPath(dPath_optimal,splineData,8);

%% linearPath
%d_linear = linearPath(d0,d1,splineData);


%% geodesicForward (TODO)
%Nsteps = 100;
%v0 = pathVelocity(d_linear,0, splineData);
%q = geodesicForward(d0,d0,Nsteps,splineData,quadData);


%% Karcher
%d_array = {d0,d1,d2};
%[d_karcher,vel,grad_norm] = karcherMeanAmpl(d_array,splineData,quadData,quadDataTensor);