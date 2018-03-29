%% curveRiemH2InnerProd
%
% Computes the Riemannian H2 inner product at the curve d and its
% derivative, i.e. it computes the derivative of the function
%   d --> G_d(v, w)
% for fixed tangent vectors v, w.
%
% Input
%   d
%       The curve. Has dimensions [N, dSpace]
%   v, w
%       Tangent vectors.
%   splineData
%       splineData describing the curve
%
% Output
%   G
%       Value of the Riemannian metric
% comp
%       The five components of the metric separately
%           [ L2, H1, H2, H1v, H1n ]
%       We have the identity
%           G = L2 + H1 + H2 + H1v + H1n
%  dG
%       Calculates the derivative of d --> G_d(v,w)
%
function [G, comp, dG] = curveRiemH2InnerProd( d, v, w, splineData )

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
Cu = quadData.Bu_S * d;
Cuu = quadData.Buu_S * d;

V = quadData.B_S * v;
Vu = quadData.Bu_S * v;
Vuu = quadData.Buu_S * v;

W = quadData.B_S * w;
Wu = quadData.Bu_S * w;
Wuu = quadData.Buu_S * w;

%% Building blocks for Riemannian metric
CuCuu = sum(Cu.*Cuu,2);
CuCuu2 = CuCuu .* CuCuu;

VW = sum(V.*W, 2);
VuWu = sum(Vu.*Wu, 2);
VuCu = sum(Vu.*Cu,2);
WuCu = sum(Wu.*Cu,2);
VuuWu = sum(Vuu.*Wu, 2);
VuWuu = sum(Vu.*Wuu, 2);
VuuWuu = sum(Vuu.*Wuu, 2);

Cspeed = sum( Cu.^2 , 2).^(1/2);
CspeedInv = 1 ./ Cspeed;
CspeedInv2 = CspeedInv .* CspeedInv;
CspeedInv3 = CspeedInv .* CspeedInv2;
CspeedInv5 = CspeedInv3 .* CspeedInv2;
CspeedInv7 = CspeedInv5 .* CspeedInv2;
CspeedInv9 = CspeedInv7 .* CspeedInv2;

%% Adjust for scale-invariance
if scaleInv
    ell = quadWeights' * Cspeed; %Length of c
    % Update coeffecients with the length weights
    a_orig = a;
    a(1) = a(1)/ell^3;
    a(2) = a(2)/ell;
    a(3) = a(3)*ell;
    a(4) = a(4)/ell;
    a(5) = a(5)/ell;
end

%% Building the metric
Ct_L2 = VW .* Cspeed;
Ct_H1 = VuWu .* CspeedInv;
Ct_H1v = VuCu .* WuCu .* CspeedInv3;
Ct_H1n = Ct_H1-Ct_H1v;
Ct_H2 = VuWu .* CuCuu.^2 .* CspeedInv7 ...
    - (VuWuu + VuuWu) .* CuCuu .* CspeedInv5 ...
    + VuuWuu .* CspeedInv3;

energyIntegrand = a(1) * Ct_L2 + a(2) * Ct_H1 + a(3) * Ct_H2 ...
    + a(4) * Ct_H1v + a(5) * Ct_H1n;

% Compute final energy
G = quadWeights' * energyIntegrand;

if nargout < 2
    return
end

comp= [quadWeights' * a(1) * Ct_L2, quadWeights' * a(2) * Ct_H1,...
       quadWeights' *a(3) * Ct_H2 ,quadWeights' *a(4) * Ct_H1v,...
       quadWeights' *a(5) * Ct_H1n]; 

%% Compute the gradient
if nargout < 3
    return
end

weight_L2 = CspeedInv .* VW .* quadWeights;
sparseW_L2 = sparse(1:noQuadSites, 1:noQuadSites, weight_L2);

weight_H1 = -CspeedInv3 .* VuWu .* quadWeights;
sparseW_H1 = sparse(1:noQuadSites, 1:noQuadSites, weight_H1);

weight_H1v_1 = CspeedInv3 .* VuCu .* quadWeights;
weight_H1v_2 = CspeedInv3 .* WuCu .* quadWeights;
weight_H1v_3 = -3 * CspeedInv5 .* VuCu .* WuCu .* quadWeights;

sparseW_H1v_1 = sparse(1:noQuadSites, 1:noQuadSites, weight_H1v_1);
sparseW_H1v_2 = sparse(1:noQuadSites, 1:noQuadSites, weight_H1v_2);
sparseW_H1v_3 = sparse(1:noQuadSites, 1:noQuadSites, weight_H1v_3);

weight_H2_1 = ( -7*CspeedInv9.*CuCuu2.*VuWu + ...
                5*CspeedInv7.*CuCuu .* (VuWuu + VuuWu) + ...
                -3*CspeedInv5.*VuuWuu ) .* quadWeights;
weight_H2_2 = ( 2*CspeedInv7.*CuCuu.*VuWu + ...
                -CspeedInv5.* (VuWuu + VuuWu) ) .* quadWeights;
sparseW_H2_1 = sparse(1:noQuadSites, 1:noQuadSites, weight_H2_1);
sparseW_H2_2 = sparse(1:noQuadSites, 1:noQuadSites, weight_H2_2);

dG_L2 = zeros([N, dSpace]);
dG_H1 = zeros([N, dSpace]);
dG_H1v = zeros([N, dSpace]);
dG_H1n = zeros([N, dSpace]);
dG_H2 = zeros([N, dSpace]);

for kk = dSpace:-1:1
    dG_L2(:,kk) = Cu(:,kk)' * sparseW_L2 * quadData.Bu_S;
    dG_H1(:,kk) = Cu(:,kk)' * sparseW_H1 * quadData.Bu_S;
    dG_H1v(:,kk) = Wu(:,kk)' * sparseW_H1v_1 * quadData.Bu_S + ...
                   Vu(:,kk)' * sparseW_H1v_2 * quadData.Bu_S + ...
                   Cu(:,kk)' * sparseW_H1v_3 * quadData.Bu_S;
    dG_H1n(:,kk) = dG_H1(:,kk) - dG_H1v(:,kk);
    dG_H2(:,kk) = Cu(:,kk)' * sparseW_H2_1 * quadData.Bu_S + ...
                    Cuu(:,kk)' * sparseW_H2_2 * quadData.Bu_S + ...
                    Cu(:,kk)' * sparseW_H2_2 * quadData.Buu_S;
end

dG = a(1) * dG_L2 + a(2) * dG_H1 + a(3) * dG_H2 + ...
        a(4) * dG_H1v + a(5) * dG_H1n;

if scaleInv
    % Derivatives of the coefficients with respect to length
    ap = zeros([5, 1]);
    if scaleInv
        ap(1) = -3 * a_orig(1) / ell^4;
        ap(2) = -a_orig(2) / ell^2;
        ap(3) = a_orig(3);
        ap(4) = -a_orig(4) / ell^2;
        ap(5) = -a_orig(5) / ell^2;
    end

    energyIntegrand2 = ap(1) * Ct_L2 + ap(2) * Ct_H1 + ap(3) * Ct_H2 ...
        + ap(4) * Ct_H1v + ap(5) * Ct_H1n;
    G2 = quadWeights' * energyIntegrand2;

    weight_ell = CspeedInv .* quadWeights;
    sparseW_ell = sparse(1:noQuadSites, 1:noQuadSites, weight_ell);

    dell = zeros([N, dSpace]);
    for kk = dSpace:-1:1
        dell(:,kk) = Cu(:,kk)' * sparseW_ell * quadData.Bu_S;
    end
    
    dG = dG + G2 * dell;
end




