%% Test derivative
splineData = constructEmptySplineData;
splineData.N = 6;
splineData.nS = 3;
splineData.Nphi = 4;
splineData.nPhi = 3;
splineData.quadDegree = [6, 4];
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);

d1 = loadDataSet('basic', splineData, '', 'curves', 'O');
ga = struct( 'phi', 0.0 * rand(splineData.Nphi, 1) , ...
             'beta', 0.3, 'v', [0.2; -0.1], 'alpha', 0.);
         
applyDiff = true;
applyTra = true;
applyRot = true;
applyShift = true;
         
[c, dc] = curveApplyGamma(d1, ga, splineData, ...
                    'applyDiff', applyDiff, 'applyTra', applyTra, ...
                    'applyRot', applyRot, 'applyShift', applyShift);

funAdimat = @(d1, phi, v, beta, alpha) ...
   curveApplyGamma(d1, struct('phi', phi, 'v', v, ...
                              'beta', beta, 'alpha', alpha), ...
                   splineData, ...
                   'applyDiff', applyDiff, 'applyTra', applyTra, ...
                   'applyRot', applyRot, 'applyShift', applyShift );
         
adopts = admOptions();
adopts.independents = [2, 3, 4, 5];

dc_AD = admDiffFD( funAdimat, 1, d1, ga.phi, ga.v, ga.beta, ga.alpha, ...
         adopts );

tol=1e-6;
% disp(['Error=', num2str(norm(dc_AD-dc))]);
assert(norm(dc_AD-dc) < tol);

