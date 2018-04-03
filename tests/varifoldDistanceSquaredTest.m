c0 = importdata('cc003.txt');
c1 = importdata('cc016.txt');

splineData = constructSplineData;
splineData.N = 40;
splineData.nS = 3;
splineData.quadDegree = [6, 4];
splineData.a = [1 2 3 4 5];
splineData.curveClosed = 1;
splineData = constructKnots(splineData);

splineData.varData = constructVarData(splineData, 'noPts', 100);
splineData = setupQuadData(splineData);

N = splineData.N;
dSpace = 2;

d0 = constructSplineApproximation(c0, splineData);
d1 = constructSplineApproximation(c1, splineData);

tau = 1e-5;

[~, dE0, dE1] = varifoldDistanceSquared(d0, d1, splineData);
%[~, dE1] = varifoldDistanceSquared(d1, d0, splineData);
% dE1 = dE0;

res0 = zeros(N, dSpace);
res1 = zeros(N, dSpace);

for ii = 1:dSpace
    for jj = 1:N
        dh = zeros(N, dSpace);
        dh(jj, ii) = 1;
    
        E0p = varifoldDistanceSquared(d0 + tau * dh, d1, splineData);
        E0m = varifoldDistanceSquared(d0 - tau * dh, d1, splineData);
        
        E1p = varifoldDistanceSquared(d0, d1 + tau * dh, splineData);
        E1m = varifoldDistanceSquared(d0, d1 - tau * dh, splineData);
        
        res01 = sum(dE0(:) .* dh(:));
        res02 = (E0p - E0m) / (2 * tau);
        
        res11 = sum(dE1(:) .* dh(:));
        res12 = (E1p - E1m) / (2 * tau);
        
        %disp(abs(res1-res2));
        res0(ii, jj) = abs(res01-res02) / max(abs(res01), abs(res02));
        res1(ii, jj) = abs(res11-res12) / max(abs(res11), abs(res12));
    end
end

disp(max(max(res0)));
disp(max(max(res1)))