function E = energyH2adimat(d1, dPathInt, d2, phi, v, beta, alpha, ...
                            splineData, optDiff, optTra, optRot, optShift)

if ~optDiff
    phi = [];
end
if ~optTra
    v = [];
end
if ~optRot
    beta = [];
end
if ~optShift
    alpha = [];
end
    
E = energyH2Diff([d1; dPathInt; d2], phi, v, beta, alpha, splineData, ...
    'optDiff', optDiff, 'optTra', optTra, 'optRot', optRot, ...
    'optShift', optShift);