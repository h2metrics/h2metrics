%% constructKnots
%
% Function constructs the knot sequences with the given parameters from
% splineData. The knot sequences are uniform in space and time, periodic in
% space and with full multiplicity at the boundary in time.
%
% interpolS is set automatically using N and nS, unless set otherwise.
%
% Input
%   splineData
%       Contains information about spline degree and number of control
%       points.
%
% Output
%   splineData
%       Same as input with knot sequences.
%
function [ splineData ] = constructKnots( splineData )

if ~isempty(splineData.N) && ~isempty(splineData.nS)
    N = splineData.N;
    nS = splineData.nS;
    % Normalize, domain of definition is [0,2*pi]
    splineData.knotsS =  [(-nS):(N+nS)]/N*2*pi; 
    splineData.innerKnotsS = splineData.knotsS(nS+1:end-nS);
    
    splineData.noInterpolS = splineData.N;
    innerKnotsS = splineData.innerKnotsS;
    if mod(splineData.nS, 2) == 0 % See deBoor for reasons for this.
        splineData.interpolS = innerKnotsS(1:end-1)' + ...
                                0.5*diff(innerKnotsS)';
    else
        splineData.interpolS = innerKnotsS(1:end-1)';
    end
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
    % Normalize, domain of definition is [0,2*pi]
    splineData.knotsPhi =  [(-nPhi):(Nphi+nPhi)]/Nphi*2*pi; 
    splineData.innerKnotsPhi = splineData.knotsPhi(nPhi+1:end-nPhi);
end

end

