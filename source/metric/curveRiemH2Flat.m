%% curveRiemH2Flat
%
% Computes the discrete flat map associated to the Riemannian metric, i.e.
% the spline p = Flat(d, v) satisfies
%    G_d(v, w) = < p, w >
% where <.,.> is the Euclidean inner product on the control points.
%
% Input
%   d
%       The curve. Has dimensions [N, dSpace]
%   v
%       Tangent vector.
%   splineData
%       splineData describing the curve
%
% Output
%   p
%       Result of the flat map applied to v at d.
%
function [p] = curveRiemH2Flat( d, v, splineData )

%% Extract parameters
a = splineData.a;
N = splineData.N;
dSpace = splineData.dSpace;
quadData = splineData.quadData;

noQuadSites = quadData.noQuadPointsS;
quadWeights = quadData.quadWeightsS;

scaleInv = 0;
if ~isempty(splineData.scaleInv)
    scaleInv = splineData.scaleInv;
end

%% Evaluate path at quadrature sites
B = quadData.B_S;
Bu = quadData.Bu_S;
Buu = quadData.Buu_S;

Cu = Bu * d;
Cuu = Buu * d;

V = quadData.B_S * v;
Vu = quadData.Bu_S * v;
Vuu = quadData.Buu_S * v;

%% Building blocks for Riemannian metric
CuCuu = sum(Cu.*Cuu,2);
CuCuu2 = CuCuu .* CuCuu;

Cspeed = sum( Cu.^2 , 2).^(1/2);
CspeedInv = 1 ./ Cspeed;
CspeedInv2 = CspeedInv .* CspeedInv;
CspeedInv3 = CspeedInv .* CspeedInv2;
CspeedInv5 = CspeedInv3 .* CspeedInv2;
CspeedInv7 = CspeedInv5 .* CspeedInv2;

VuCu = sum(Vu.*Cu,2);

%% Building the momentum
weight_L2 = Cspeed .* quadWeights;
sparseW_L2 = sparse(1:noQuadSites, 1:noQuadSites, weight_L2);

weight_H1 = CspeedInv .* quadWeights;
sparseW_H1 = sparse(1:noQuadSites, 1:noQuadSites, weight_H1);

weight_H1v = VuCu .* CspeedInv3 .* quadWeights;
sparseW_H1v = sparse(1:noQuadSites, 1:noQuadSites, weight_H1v);

weight_H2_1 = CspeedInv7.*CuCuu2 .* quadWeights;
weight_H2_2 = -CspeedInv5.*CuCuu .* quadWeights;
weight_H2_3 = CspeedInv3 .* quadWeights;
sparseW_H2_1 = sparse(1:noQuadSites, 1:noQuadSites, weight_H2_1);
sparseW_H2_2 = sparse(1:noQuadSites, 1:noQuadSites, weight_H2_2);
sparseW_H2_3 = sparse(1:noQuadSites, 1:noQuadSites, weight_H2_3);

p_L2 = zeros([N, dSpace]);
p_H1 = zeros([N, dSpace]);
p_H1v = zeros([N, dSpace]);
p_H1n = zeros([N, dSpace]);
p_H2 = zeros([N, dSpace]);

for kk = dSpace:-1:1
    p_L2(:,kk) = V(:,kk)' * sparseW_L2 * B;
    p_H1(:,kk) = Vu(:,kk)' * sparseW_H1 * Bu;
        
    p_H2(:,kk) = Vu(:,kk)' * sparseW_H2_1 * Bu + ...
                 Vuu(:,kk)' * sparseW_H2_2 * Bu + ...
                 Vu(:,kk)' * sparseW_H2_2 * Buu + ...
                 Vuu(:,kk)' * sparseW_H2_3 * Buu;
                     
    p_H1v(:,kk) = Cu(:,kk)' * sparseW_H1v * Bu;
    p_H1n(:,kk) = p_H1(:,kk) - p_H1v(:,kk);                     
end

%% Scale-invariant metrics
if scaleInv
    ell = quadWeights' * Cspeed; % Length of c
    
    % Update coeffecients with the length weights
    a(1) = a(1)/ell^3;
    a(2) = a(2)/ell;
    a(3) = a(3)*ell;
    a(4) = a(4)/ell;
    a(5) = a(5)/ell;
end

%% Put it together
p = a(1)*p_L2 + a(2)*p_H1 + a(3)*p_H2 + a(4)*p_H1v + a(5)*p_H1n;

end