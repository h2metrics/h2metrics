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
%   splineData, quadData
%       Contains interpolation sites and collocation matrices.
%
% Output
%   c
%       Spline curve approximating d o phi.
function c = curveComposeDiff(d, phi, splineData, quadData)

nS = splineData.nS;
knotsS = splineData.knotsS;

B_interpolS = quadData.B_interpolS;
B_interpolPhi = quadData.B_interpolPhi;
     
phiPts = B_interpolPhi * phi;
phiPts = phiPts + splineData.interpolS;
phiPts = mod(phiPts, 2*pi);

% fastBsplineEval
% dSpace = splineData.dSpace;
% dNonper = [ d; d(1:nS,:) ];
% for ii = dSpace:-1:1
%     cPts(:,ii) = fastBSplineEval(knotsS, dNonper(:,ii), nS, phiPts);
% end

% For compatibility purposes, use deBoor.
cPts = deBoor( knotsS, nS, d, phiPts, 1, 'periodic', true );

c = B_interpolS \ cPts;

% When we wanted to use first derivatives for interpolation...
% [ knots_dx,weights_dx,order_dx] = fastBSplineDer( knotsS,dNonper(:,1),nS );
% [ knots_dx,weights_dy,order_dx] = fastBSplineDer( knotsS,dNonper(:,2),nS );
% dcPts(:,1) = fastBSplineEval(knots_dx, weights_dx, nS-1, phiPts);
% dcPts(:,2) = fastBSplineEval(knots_dx, weights_dy, nS-1, phiPts);

% knots_dx_alt = knotsS(2:end-1);
% weights_dx_alt = weights_dx(2:end-1);
% weights_dy_alt = weights_dy(2:end-1);
% dcPts_alt(:,1) = fastBSplineEval(knots_dx_alt, weights_dx_alt, nS-1, phiPts);
% dcPts_alt(:,2) = fastBSplineEval(knots_dx_alt, weights_dy_alt, nS-1, phiPts);

% [ knots_dx,weights_dx,order_dx] = fastBSplineDer( knotsS,dNonper(:,1),nS )

end
