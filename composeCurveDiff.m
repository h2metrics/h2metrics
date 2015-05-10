%% composeCurveDiff
%
% Function computes the composition of a periodic curve and a
% diffeomorphism. It uses the fastBSpline library to speed up spline
% evaluation.
%
% Points that are used to fit the control points of the composition are
% given by splineData.interpolS.
%
% Input
%   d
%       Curve
%   phi
%       Diffeomorphism
%   splineData, quadData
%       Contains interpolation sites and collocation matrices.
%
% Output
%   c
%       Spline curve approximating d o phi.
function c = composeCurveDiff(d, phi, splineData, quadData)

dSpace = splineData.dSpace;
nS = splineData.nS;
knotsS = splineData.knotsS;

B_interpolS = quadData.B_interpolS;
B_interpolPhi = quadData.B_interpolPhi;
     
phiPts = B_interpolPhi * phi;
phiPts = phiPts + splineData.interpolS;

dNonper = [ d; d(1:nS,:) ];

for ii = dSpace:-1:1
    cPts(:,ii) = fastBSplineEval(knotsS, dNonper(:,ii), nS, phiPts);
end

c = B_interpolS \ cPts;

end