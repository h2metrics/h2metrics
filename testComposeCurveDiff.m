% Script for test testing all functionality
%% constructEmptySplineData, and setup parameters
splineData = constructEmptySplinedata;
splineData.N = 100; %no. control points, must be bigger than n+1
splineData.Nphi = 100; %No. control points for diffeomorphisms
splineData.nS = 5; %spacial degree
splineData.nPhi = 3; %diffemorphism degree
splineData.dSpace = 2;
splineData.noInterpolS = 2 * splineData.N; % For composition

splineData.Nt = 7 + 2; %Number of time control points
splineData.nT = 1; %time degree
splineData.quadDegree = [8,4]; %Quadrature precission

%% constructKnots and quadData
splineData = constructKnots(splineData);
[quadData, quadDataTensor] = setupQuadData(splineData);

%% constructSplineApproximation
f0 = @(t) 1/100*[17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - 2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + 13*sin(3*t) - sin(4*t) - 5*sin(5*t),...
 -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + 13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + 19*sin(3*t) - 8*sin(4*t)];
% f1 = @(t) 1/100*[150*cos(-t - 3*pi/4), 150*sin(-t - 3*pi/4)] ;

g = @(t) 0.3 * sin(3*t);
f2 = @(t) f0(t + g(t));

d0 = constructSplineApproximation(f0, splineData);
d2 = constructSplineApproximation(f2, splineData);

interpolS = splineData.interpolS;
B_interpolPhi = quadData.B_interpolPhi;
phi = B_interpolPhi \ g(interpolS);

%% Test composition
c2 = composeCurveDiff(d0, phi, splineData, quadData);

B_interpolS = quadData.B_interpolS;
y1 = B_interpolS * d2;
y2 = B_interpolS * c2;
disp(max(max(abs(y1-y2))));