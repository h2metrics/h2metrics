%% Check that it coincides with curveApplyGamma
% optDiff=true; optTra=true; optRot=true; optShift=true

splineData = constructEmptySplineData;
splineData.N = 55; % Problem for N=55 and higher,
                   % works if replaced by 54 or lower
splineData.nS = 4;
splineData.Nt = 15 + 2;
splineData.nT = 2;
splineData.Nphi = 10;
splineData.nPhi = 4;
splineData.quadDegree = [6, 4];
splineData = constructKnots(splineData);
[quadData, quadDataTensor, splineData] = setupQuadData(splineData);
splineData.a = [1 0 1e-4];

d1 = loadDataSet('basic', splineData, '', 'curves', 'U');
d2 = loadDataSet('basic', splineData, '', 'curves', 'T');
ga = struct( 'phi', 0.01 * rand(splineData.Nphi, 1) , ...
             'beta', 0.3, 'v', [0.2; -0.1], 'alpha', [0.1]);
d2Ga = curveApplyGamma(d2, ga, splineData, quadData);
dPathGa = linearPath(d1, d2Ga, splineData);

% Result 1 - use energyH2Diff as in geodesicBvpParam
res1 = energyH2Diff( ...
    [d1; dPathGa(splineData.N+1:end-splineData.N, :); d2], ...
    ga.phi, ga.v, ga.beta, ga.alpha, ...
    splineData, quadData, quadDataTensor, ...
    'optDiff', true, 'optTra', true, 'optRot', true, ...
    'optShift', true );

% Result 2 - apply gamma manually
res2 = energyH2Diff( dPathGa, [], [], [], [], ...
    splineData, quadData, quadDataTensor, ...
    'optDiff', false, 'optTra', false, 'optRot', false, ...
    'optShift', false );

% Result 3 - pathRiemH2Energy
res3 = pathRiemH2Energy(dPathGa, splineData, quadData, quadDataTensor);

tol = 1e-12;
assert(abs(res1-res2) < tol);
assert(abs(res2-res3) < tol);
assert(abs(res3-res1) < tol);

%% Test gradient, no groups
% optDiff=false; optTra=false; optRot=false; optShift=false

splineData = constructEmptySplineData;
splineData.N = 10; % Problem for N=55 and higher,
                   % works if replaced by 54 or lower
splineData.nS = 4;
splineData.Nt = 10 + 2;
splineData.nT = 2;
splineData.Nphi = 7;
splineData.nPhi = 4;
splineData.quadDegree = [6, 4];
splineData = constructKnots(splineData);
[quadData, quadDataTensor, splineData] = setupQuadData(splineData);
splineData.a = [1 0 1e-4];

d1 = loadDataSet('basic', splineData, '', 'curves', 'H');
d2 = loadDataSet('basic', splineData, '', 'curves', 'O');
dPath = linearPath(d1, d2, splineData);
dPathInt = dPath(splineData.N+1:end-splineData.N,:);

optDiff = false;
optTra = false;
optRot = false;
optShift = false;

[E, dE, H] = energyH2Diff(dPath, [], [], [], [], ...
    splineData, quadData, quadDataTensor, ...
    'optDiff', optDiff, 'optTra', optTra, 'optRot', optRot, ...
    'optShift', optShift);

adopts = admOptions();
adopts.functionResults = {E}; % for admDiffRev
adopts.independents = [2, 4, 5, 6, 7];

dE_AD = admDiffFD( @energyH2adimat, 1, d1, dPathInt, d2, ...
         zeros(1, 1), [0; 0], 0, 0, ...
         splineData, optDiff, optTra, optRot, optShift, adopts );
dE_AD = dE_AD';

tol=1e-6;
assert(norm(dE_AD-dE) < tol);

%% Test gradient, all groups
% optDiff=true; optTra=true; optRot=true; optShift=true

splineData = constructEmptySplineData;
splineData.N = 33; % Problem for N=55 and higher,
                   % works if replaced by 54 or lower
splineData.nS = 4;
splineData.Nt = 15 + 2;
splineData.nT = 2;
splineData.Nphi = 12;
splineData.nPhi = 4;
splineData.quadDegree = [6, 4];
splineData = constructKnots(splineData);
[quadData, quadDataTensor, splineData] = setupQuadData(splineData);
splineData.a = [1 0 1e-4];

d1 = loadDataSet('basic', splineData, '', 'curves', 'H');
d2 = loadDataSet('basic', splineData, '', 'curves', 'O');
ga = struct( 'phi', 0.0 * rand(splineData.Nphi, 1) , ...
             'beta', 0.3, 'v', [0.2; -0.1], 'alpha', [0.1]);
d2Ga = curveApplyGamma(d2, ga, splineData, quadData);
dPathGa = linearPath(d1, d2Ga, splineData);
dPathGaInt = dPathGa(splineData.N+1:end-splineData.N,:);
dPathGa2 = [d1; dPathGaInt; d2];

optDiff = true;
optTra = true;
optRot = true;
optShift = true;

[E, dE, H] = energyH2Diff(dPathGa2, ga.phi, ga.v, ga.beta, ga.alpha, ...
    splineData, quadData, quadDataTensor, ...
    'optDiff', optDiff, 'optTra', optTra, 'optRot', optRot, ...
    'optShift', optShift);

adopts = admOptions();
adopts.functionResults = {E}; % for admDiffRev
adopts.independents = [2, 4, 5, 6, 7];

dE_AD = admDiffFD( @energyH2adimat, 1, d1, dPathGaInt, d2, ...
         ga.phi, ga.v, ga.beta, ga.alpha, ...
         splineData, optDiff, optTra, optRot, optShift, adopts );
dE_AD = dE_AD';

tol=1e-3;
%disp(['N=', num2str(jj), ' Error=', num2str(norm(dE_AD-dE))]);
% norm(dE_AD-dE)
assert(norm(dE_AD-dE) < tol);


%% Test gradient, all group combinations
% optDiff=true; optTra=true; optRot=true; optShift=true

optDiffList =  [false, false, false, false, false, false, false, false, ...
                true,  true,  true,  true,  true,  true,  true,  true];
optTraList =   [false, false, false, false, true,  true,  true,  true, ...
                false, false, false, false, true,  true,  true,  true];
optRotList =   [false, false, true,  true,  false, false, true,  true, ...
                false, false, true,  true,  false, false, true,  true];
optShiftList = [false, true,  false, true,  false, true,  false, true, ...
                false, true,  false, true,  false, true,  false, true];
            
for kk = 1:16
    splineData = constructEmptySplineData;
    splineData.N = 13;
    splineData.nS = 3;
    splineData.Nt = 15 + 2;
    splineData.nT = 2;
    splineData.Nphi = 4;
    splineData.nPhi = 3;
    splineData.quadDegree = [6, 4];
    splineData = constructKnots(splineData);
    [quadData, quadDataTensor, splineData] = setupQuadData(splineData);
    splineData.a = [1 0 1e-4];

    d1 = loadDataSet('basic', splineData, '', 'curves', 'H');
    d2 = loadDataSet('basic', splineData, '', 'curves', 'O');
    ga = struct( 'phi', 0.0 * rand(splineData.Nphi, 1) , ...
                 'beta', 0.3, 'v', [0.2; -0.1], 'alpha', [0.1]);
    d2Ga = curveApplyGamma(d2, ga, splineData, quadData);
    dPathGa = linearPath(d1, d2Ga, splineData);
    dPathGaInt = dPathGa(splineData.N+1:end-splineData.N,:);
    dPathGa2 = [d1; dPathGaInt; d2];

    optDiff = optDiffList(kk);
    optTra = optTraList(kk);
    optRot = optRotList(kk);
    optShift = optShiftList(kk);
    
    if ~optDiff
        ga.phi = 0;
    end

    [E, dE, H] = energyH2Diff(dPathGa2, ga.phi, ga.v, ga.beta, ga.alpha, ...
        splineData, quadData, quadDataTensor, ...
        'optDiff', optDiff, 'optTra', optTra, 'optRot', optRot, ...
        'optShift', optShift);

    adopts = admOptions();
    %adopts.functionResults = {E}; % for admDiffRev
    adopts.independents = [2, 4, 5, 6, 7];

    dE_AD = admDiffFD( @energyH2adimat, 1, d1, dPathGaInt, d2, ...
             ga.phi, ga.v, ga.beta, ga.alpha, ...
             splineData, optDiff, optTra, optRot, optShift, adopts );
    dE_AD = dE_AD';

    tol=1e-5;
    disp(['D=', num2str(optDiff), ' T=', num2str(optTra), ...
          ' R=', num2str(optRot), ' S=', num2str(optShift), ...
          ' Error=', num2str(norm(dE_AD-dE))]);
    assert(norm(dE_AD-dE) < tol);
end
