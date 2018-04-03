%% Test 1 - H2 metric
f0 = @(t) 1/100*[ 17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - ...
                  2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + ...
                  13*sin(3*t) - sin(4*t) - 5*sin(5*t), ...
                  -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
                  13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + ...
                  19*sin(3*t) - 8*sin(4*t)];
v0 = @(t) [3*cos(t) - 5*cos(2*t) + sin(4*t), ...
          2*cos(2*t) - 3*sin(4*t)];
w0 = @(t) [2*cos(t) - 6*cos(3*t) + sin(2*t), ...
          cos(2*t) + 2*sin(3*t)];
tol = 1e-12;

splineData = constructSplineData;
splineData.N = 21;
splineData.nS = 3;
splineData.quadDegree = [6, 4];
splineData.a = [2, 1, 3, 0.5, 7];
splineData.curveClosed = 1;
splineData.scaleInv = 0;
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);

d = constructSplineApproximation(f0, splineData);
v = constructSplineApproximation(v0, splineData);
w = constructSplineApproximation(w0, splineData);

[G, ~, dG] = curveRiemH2InnerProd(d, v, w, splineData);
dG1v = sum(dG(:) .* w(:));

tau = 1e-7;
dp = d + tau * w;
dm = d - tau * w;
dG2v = (curveRiemH2InnerProd(dp, v, w, splineData) - ...
        curveRiemH2InnerProd(dm, v, w, splineData)) / (2*tau);



