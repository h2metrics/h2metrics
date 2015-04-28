function [ splineData ] = constructKnots( splineData )
% Given parameters for the the BVP problem, compute knot vectors
% Knots are assumed to be uniform.

splineData.knotsS =  [(-splineData.nS):(splineData.N+splineData.nS)]/splineData.N*2*pi; %normalize, domain of definition is [0,2*pi]x[0,2*pi]
splineData.knotsT = [ zeros(1,splineData.nT), linspace(0,1,splineData.Nt- splineData.nT + 1), ones(1,splineData.nT)];
splineData.knotsPhi =  [(-splineData.nPhi):(splineData.Nphi+splineData.nPhi)]/splineData.Nphi*2*pi; %normalize, domain of definition is [0,2*pi]x[0,2*pi]

splineData.innerKnotsS = splineData.knotsS(splineData.nS+1:end-splineData.nS);
splineData.innerKnotsT = splineData.knotsT(splineData.nT+1:end-splineData.nT);
splineData.innerKnotsPhi = splineData.knotsPhi(splineData.nPhi+1:end-splineData.nPhi);

end

