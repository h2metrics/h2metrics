%% Test 1 - function Handle
f = @(t) 1/100*[ 17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - ...
                 2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + ...
                 13*sin(3*t) - sin(4*t) - 5*sin(5*t), ...
                 -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
                 13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + ...
                 19*sin(3*t) - 8*sin(4*t)];
              
tol = 1e-14;

splineData = constructEmptySplineData;
splineData.N = 30; %no. control points, must be bigger than n+1
splineData.nS = 4; %spacial degree
splineData.curveClosed = 1;
splineData = constructKnots(splineData);

d = constructSplineApproximation(f, splineData);

th = splineData.interpolS;
pts1 = evalCurve(th, d, splineData);
pts2 = f(th);

assert(norm(pts1-pts2) < tol);

%% Test 2 - set of points
numPts = 100;
t = linspace(0, 2*pi, numPts+1)';
t = t(1:end-1);

fun = @(t) 1/100*[ 17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - ...
                   2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + ...
                   13*sin(3*t) - sin(4*t) - 5*sin(5*t), ...
                   -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
                   13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + ...
                   19*sin(3*t) - 8*sin(4*t)];
               
tol = 1e-14;

splineData = constructEmptySplineData;
splineData.N = 30; %no. control points, must be bigger than n+1
splineData.nS = 4; %spacial degree
splineData = constructKnots(splineData);

f = fun(t);
d = constructSplineApproximation(f, splineData);

pts1 = evalCurve(t, d, splineData);
pts2 = f;

assert(norm(pts1-pts2) < 1e-3);

%% Test 3 - function handle, open curve
f = @(t) 1/100*[ t + t.^2, ...
                 -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
                 13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + ...
                 19*sin(3*t) - 8*sin(4*t)];
              
tol = 1e-14;

splineData = constructEmptySplineData;
splineData.N = 30;
splineData.nS = 4; 
splineData.curveClosed = 0;
splineData = constructKnots(splineData);

d = constructSplineApproximation(f, splineData);

th = splineData.interpolS;
pts1 = evalCurve(th, d, splineData);
pts2 = f(th);

assert(norm(pts1-pts2) < tol);

%% Test 4 - set of points, open curve
numPts = 100;
t = linspace(0, 2*pi, numPts)';

fun = @(t) 1/100*[ t + t.^2, ...
                   -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
                   13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + ...
                   19*sin(3*t) - 8*sin(4*t)];
               
tol = 1e-14;

splineData = constructEmptySplineData;
splineData.N = 40;
splineData.nS = 3; 
splineData.curveClosed = 0;
splineData = constructKnots(splineData);

f = fun(t);
d = constructSplineApproximation(f, splineData);

pts1 = evalCurve(t, d, splineData);
pts2 = f;

assert(norm(pts1-pts2) < 1e-3);
