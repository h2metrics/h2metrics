%% Test 1 - H2 metric
f0 = @(t) 1/100*[ 17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - ...
                  2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + ...
                  13*sin(3*t) - sin(4*t) - 5*sin(5*t), ...
                  -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
                  13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + ...
                  19*sin(3*t) - 8*sin(4*t)];
f1 = @(t) 1/100*[ 17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - ...
                  2*cos(4*t) - 24*sin(t) - 86*sin(2*t) + ...
                  13*sin(3*t), ...
                  -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
                  13*cos(4*t) - 11*cos(5*t) + 79*sin(t)];

tol = 1e-3;

splineData = constructSplineData;
splineData.N = 21; %no. control points, must be bigger than n+1
splineData.nS = 3; %spacial degree
splineData.quadDegree = [6, 4];
splineData.a = [1 2 3 4 5 ];
splineData.curveClosed = 1;
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);

c0 = constructSplineApproximation(f0, splineData);
c1 = constructSplineApproximation(f1, splineData);
dPath = linearPath(c0, c1, splineData);

Ln = pathRiemH2Length2(dPath, splineData);
En = pathRiemH2Energy(dPath, splineData);

G1 = Ln^2;
G2 = En;

assert(abs(G1-G2) / max(abs(G1), abs(G2)) < tol);


