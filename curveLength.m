%% curveLength
%
% Calculates the length of the curve. The formula is
%   ell(d) = \int_{S^1} |d'| d\theta
%
% Input
%   d
%       Control points of the curve
%   splineData
%       General information about the splines used. (Not used currently.)
%   quadData
%       Precomputed spline collocation matrices at quadrature points.
%
% Output
%   ell
%       Length of the curve
function ell = curveLength(d, splineData, quadData)

cQuad_u = quadData.Bu_S * d;
cSpeed = sum( cQuad_u.^2 , 2).^(1/2);
ell = sum(cSpeed .* quadData.quadWeightsS);

end