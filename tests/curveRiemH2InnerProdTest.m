%% Test1 - AdiMat
f0 = @(t) 1/100*[ 17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - ...
                  2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + ...
                  13*sin(3*t) - sin(4*t) - 5*sin(5*t), ...
                  -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
                  13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + ...
                  19*sin(3*t) - 8*sin(4*t)];
u0 = @(t) [3*cos(t) - 5*cos(2*t) + sin(4*t), ...
          2*cos(2*t) - 3*sin(4*t)];
v0 = @(t) [2*cos(t) - 6*cos(3*t) + sin(2*t), ...
           cos(2*t) + 2*sin(3*t)];
tol = 1e-10;

splineData = constructEmptySplineData;
splineData.N = 20; %no. control points, must be bigger than n+1
splineData.nS = 3; %spacial degree
splineData.quadDegree = [6, 4];
splineData.a = [1, 1, 1];
splineData = constructKnots(splineData);
quadData = setupQuadData(splineData);

d = constructSplineApproximation(f0, splineData);
u = constructSplineApproximation(u0, splineData);
v = constructSplineApproximation(v0, splineData);

% Now run adiMat
G = curveRiemH2InnerProd(d, u, v, splineData, quadData);

adopts = admOptions();
adopts.functionResults = {G}; % for admDiffRev
adopts.independents = [1];

J_AD = admDiffRev( @curveRiemH2InnerProd, 1, d, u, v, ...
                   splineData, quadData, adopts );
J_AD = reshape(J_AD, [splineData.N, splineData.dSpace]);

[~, Gder] = curveRiemH2InnerProd(d, u, v, splineData, quadData);

assert(norm(Gder-J_AD) < tol);
