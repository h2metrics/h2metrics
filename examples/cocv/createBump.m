function [ d ] = createBump( centre, width, height, splineData)

%% Create bump function
sdBump = constructSplineData;
sdBump.curveClosed = 0;
sdBump.N = 13;
sdBump.nS = 3;
sdBump = constructKnots(sdBump);
sdBump = setupQuadData(sdBump);

c0 = pi;
w0 = 2*pi / 2.5; % Width of original bump

% Spline of original bump
dBump = constructSplineApproximation(@(t) [t, 0.*t], sdBump);
dBump(7,2) = 1;

h0 = evalCurve(pi, dBump, sdBump);
h0 = h0(2);

f1 = @(x) x;
f2 = @(x) evalCurve( w0*(x-centre)/width+c0, dBump(:,2), sdBump) / h0 * height;
d(:,1) = constructSplineApproximation(f1, splineData);
d(:,2) = constructSplineApproximation(f2, splineData);




end

