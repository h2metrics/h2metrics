function [ splineData ] = constructKnots( splineData )
% Given parameters for the the BVP problem, compute knot vectors
% Knots are assumed to be uniform.

if ~isempty(splineData.N) && ~isempty(splineData.nS)
    N = splineData.N;
    nS = splineData.nS;
    % Normalize, domain of definition is [0,2*pi]
    splineData.knotsS =  [(-nS):(N+nS)]/N*2*pi; 
    splineData.innerKnotsS = splineData.knotsS(nS+1:end-nS);
end

if ~isempty(splineData.Nt) && ~isempty(splineData.nT)
    Nt = splineData.Nt;
    nT = splineData.nT;
    splineData.knotsT = [ zeros(1,nT), linspace(0,1,Nt-nT+1), ...
                          ones(1,nT) ];
    splineData.innerKnotsT = splineData.knotsT(nT+1:end-nT);
end

if ~isempty(splineData.Nphi) && ~isempty(splineData.nPhi)
    Nphi = splineData.Nphi;
    nPhi = splineData.nPhi;
    %normalize, domain of definition is [0,2*pi]
    splineData.knotsPhi =  [(-nPhi):(Nphi+nPhi)]/Nphi*2*pi; 
    splineData.innerKnotsPhi = splineData.knotsPhi(nPhi+1:end-nPhi);
end

if ~isempty(splineData.noInterpolS)
    splineData.interpolS = linspace( 0, 2*pi, splineData.noInterpolS+1)';
    % Last point correponds to first point
    splineData.interpolS = splineData.interpolS(1:end-1); 
end

end

