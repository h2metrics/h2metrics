%% curveApplyShift
%
% Applies a shift by -alpha to the curve. The formula is
%   c(th) = d( th - alpha )
%
% Input
%   d
%       Control points of the curve
%   alpha
%       Shift to be applied
%   splineData
%       General information about the splines used.%
% Output
%   c
%       The transformed curve
function c = curveApplyShift(d, alpha, splineData)

nS = splineData.nS;
knotsS = splineData.knotsS;

interpolS = splineData.interpolS;
phiPts = mod(interpolS - alpha, 2*pi);
cPts = deBoor( knotsS, nS, d, phiPts, 1, 'periodic', true );

c = splineData.quadData.B_interpolS \ cPts;

end
