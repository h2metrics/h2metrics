c0 = importdata('cc003.txt');
c1 = importdata('cc016.txt');
c2 = importdata('cc073.txt');

splineData = constructSplineData;
splineData.N = 40;
splineData.nS = 3;
splineData.Nt = 5;
splineData.nT = 2;
splineData.quadDegree = [6, 4];
splineData.a = [1 2 3 4 5];
splineData.curveClosed = 1;
splineData.scaleInv = 0;
splineData = constructKnots(splineData);

splineData.varData = constructVarData(splineData, 'noPts', 100, ...
                                      'kernelSizeGeom', 100);
splineData = setupQuadData(splineData);

N = splineData.N;
Nt = splineData.Nt;
dSpace = 2;

d0 = constructSplineApproximation(c0, splineData);
d1 = constructSplineApproximation(c1, splineData);
d2 = constructSplineApproximation(c2, splineData);

dPath = linearPath(d0, d1, splineData);
rho = 1.1;
beta = 0.01;
v = [0.1; 0.2];

coeffs = [ reshape(dPath(N+1:end, :), [], 1); ...
              rho; beta; v ];
          
pars = struct();
pars.splineData = splineData;
pars.optRot = 1;
pars.optTra = 1;
pars.optScal = 1;
pars.checkTurningNumber = 0;
pars.d0 = d0;
pars.dEnd = d2;
pars.lambda = 0.1;

tau = 1e-4;

[~, dE1] = energyH2VarifoldQuadraticPenalty(coeffs, pars, splineData);

res = zeros(size(coeffs));
for jj = 1:length(coeffs)
    dcoeffs = zeros(size(coeffs));
    dcoeffs(jj) = 1;
    
    Ep = energyH2VarifoldQuadraticPenalty(coeffs + tau * dcoeffs, ...
                                          pars, splineData);
    Em = energyH2VarifoldQuadraticPenalty(coeffs - tau * dcoeffs, ...
                                          pars, splineData);
    res1 = sum(dE1 .* dcoeffs);
    res2 = (Ep - Em) / (2 * tau);
    %disp(abs(res1-res2));
    res(jj) = abs(res1-res2) / max(abs(res1), abs(res2));
end

disp(max(res));