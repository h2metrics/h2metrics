%% curveCenter
%
% Calculates the center of mass of the curve and centers it, such that the
% center of mass is moved to the origin. The formula is
%   center(d) = length(d)^(-1) * \int_{S^1} d |d'| d\theta
% and
%   c = d - center(d)
%
% Input
%   d
%       Control points of the curve
%   splineData
%       General information about the splines used.
%   quadData
%       Precomputed spline collocation matrices at quadrature points.
%
% Output
%   c
%       The transformed curve
%   center
%       Center of mass of d
function [ c, center ] = curveCenter( d, splineData, quadData )

N = splineData.N;
dSpace = splineData.dSpace;

cQuad = quadData.B_S * d;
cQuad_u = quadData.Bu_S * d;
cSpeed = sum( cQuad_u.^2 , 2).^(1/2);
cLength = sum(cSpeed .* quadData.quadWeightsS);

for jj = dSpace:-1:1
    cCenter(jj, 1) = sum(cQuad(:,jj) .* cSpeed ...
        .* quadData.quadWeightsS) / cLength;
end

c = d - ones([N, 1]) * cCenter';
center = cCenter;

end

