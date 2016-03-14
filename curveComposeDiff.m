%% curveComposeDiff
%
% Function computes the composition of a curve and a diffeomorphism.
%
% Points that are used to fit the control points of the composition are
% given by splineData.interpolS.
%
% Input
%   d
%       Curve
%   phi
%       Diffeomorphism
%   splineData
%       Contains interpolation sites and collocation matrices.
%
% Output
%   c
%       Spline curve approximating d o phi.
function c = curveComposeDiff(d, phi, splineData)

nS = splineData.nS;
knotsS = splineData.knotsS;
quadData = splineData.quadData;

B_interpolS = quadData.B_interpolS;
B_interpolPhi = quadData.B_interpolPhi;
     
phiPts = B_interpolPhi * phi;
phiPts = phiPts + splineData.interpolS;
phiPts = mod(phiPts, 2*pi);

cPts = deBoor( knotsS, nS, d, phiPts, 1, 'periodic', true );

c = B_interpolS \ cPts;

end
