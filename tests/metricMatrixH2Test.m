%% Test 1 - H2 metric
f0 = @(t) 1/100*[ 17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - ...
                  2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + ...
                  13*sin(3*t) - sin(4*t) - 5*sin(5*t), ...
                  -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
                  13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + ...
                  19*sin(3*t) - 8*sin(4*t)];
% u0 = @(t) [3*cos(t) - 5*cos(2*t) + sin(4*t), ...
%           2*cos(2*t) - 3*sin(4*t)];
% v0 = @(t) [2*cos(t) - 6*cos(3*t) + sin(2*t), ...
%           cos(2*t) +2*sin(3*t)];
tol = 1e-10;

splineData = constructEmptySplineData;
splineData.N = 21; %no. control points, must be bigger than n+1
splineData.nS = 3; %spacial degree
splineData.quadDegree = [6, 4];
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);

d = constructSplineApproximation(f0, splineData);

matG = metricMatrixH2(d, splineData);

% Now compute by hand
N = splineData.N;
dSpace = splineData.dSpace;

matGalt = [];
for jj = N*dSpace:-1:1
    for kk = N*dSpace:-1:jj
        U = zeros([N*dSpace, 1]);
        V = zeros([N*dSpace, 1]);
        U(jj) = 1;
        V(kk) = 1;
        
        u = reshape(U, [N, dSpace]);
        v = reshape(V, [N, dSpace]);
        
        matGalt(jj, kk) = curveRiemH2InnerProd(d, u, v, splineData);
        matGalt(kk, jj) = matGalt(jj, kk);
    end
end

assert(norm(matG - matGalt) < tol);


%% Test 2 - H1 metric
f0 = @(t) 1/100*[ 17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - ...
                  2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + ...
                  13*sin(3*t) - sin(4*t) - 5*sin(5*t), ...
                  -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
                  13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + ...
                  19*sin(3*t) - 8*sin(4*t)];
tol = 1e-10;

splineData = constructEmptySplineData;
splineData.N = 21; %no. control points, must be bigger than n+1
splineData.nS = 3; %spacial degree
splineData.quadDegree = [6, 4];
splineData.a = [1, 1, 0];
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);

d = constructSplineApproximation(f0, splineData);

matG = metricMatrixH2(d, splineData);

% Now compute by hand
N = splineData.N;
dSpace = splineData.dSpace;

matGalt = [];
for jj = N*dSpace:-1:1
    for kk = N*dSpace:-1:jj
        U = zeros([N*dSpace, 1]);
        V = zeros([N*dSpace, 1]);
        U(jj) = 1;
        V(kk) = 1;
        
        u = reshape(U, [N, dSpace]);
        v = reshape(V, [N, dSpace]);
        
        matGalt(jj, kk) = curveRiemH2InnerProd(d, u, v, splineData);
        matGalt(kk, jj) = matGalt(jj, kk);
    end
end

assert(norm(matG - matGalt) < tol);
