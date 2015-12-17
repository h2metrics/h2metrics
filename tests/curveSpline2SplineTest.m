%% Test 1 - spline to itself
f0 = @(t) 1/100*[ 17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - ...
                  2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + ...
                  13*sin(3*t) - sin(4*t) - 5*sin(5*t), ...
                  -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
                  13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + ...
                  19*sin(3*t) - 8*sin(4*t)];
tol = 1e-14;

splineData = constructEmptySplineData;
splineData.N = 10; %no. control points, must be bigger than n+1
splineData.nS = 1; %spacial degree
splineData = constructKnots(splineData);

d0 = constructSplineApproximation(f0, splineData);
d1 = curveSpline2Spline(d0, splineData, splineData);

assert(norm(d0-d1) < tol);

%% Test 2 - spline to higher order, evaluate at some points
f0 = @(t) 1/100*[ 17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - ...
                  2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + ...
                  13*sin(3*t) - sin(4*t) - 5*sin(5*t), ...
                  -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
                  13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + ...
                  19*sin(3*t) - 8*sin(4*t)];
tol = 1e-14;

splineData = constructEmptySplineData;
splineData.N = 11; %no. control points, must be bigger than n+1
splineData.nS = 2; %spacial degree
splineData = constructKnots(splineData);

splineDataNew = constructEmptySplineData;
splineDataNew.N = 500; %no. control points, must be bigger than n+1
splineDataNew.nS = 5; %spacial degree
splineDataNew = constructKnots(splineDataNew);

d1 = constructSplineApproximation(f0, splineData);
d2 = curveSpline2Spline(d1, splineData, splineDataNew);

noPts = 100;
x = linspace(0, 2*pi, noPts);
y1 = deBoor(splineData.knotsS, splineData.nS, d1, x, 1, 'periodic', true);
y2 = deBoor(splineDataNew.knotsS, splineDataNew.nS, d2, x, 1, ...
            'periodic', true);

assert(max(max(abs(y1-y2))) < 1e-3);