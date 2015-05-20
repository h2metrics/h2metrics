function G = curveRiemH2InnerProd( d, v, w, splineData, quadData )

%% Extract parameters
dSpace = splineData.dSpace;
N = splineData.N;

%% Evaluate path at quadrature sites
Cu = quadData.Bu_S * d;
Cuu = quadData.Buu_S * d;

V = quadData.B_S * v;
Vu = quadData.Bu_S * v;
Vuu = quadData.Buu_S * v;

W = quadData.B_S * w;
Wu = quadData.Bu_S * w;
Wuu = quadData.Buu_S * w;

%% Remaining stuff (should be identical to energyH2)
Cspeed = Cu(:,1) .* Cu(:,1);
for ii = 2:dSpace
    Cspeed = Cspeed + Cu(:,ii) .* Cu(:,ii);
end
Cspeed = sqrt(Cspeed);
CspeedInv = 1 ./ Cspeed;

% L2 and H1 terms
Ct_L2 = V(:,1) .* W(:,1);
Ct_H1 = Vu(:,1) .* Wu(:,1);
for ii = 2:dSpace
    Ct_L2 = Ct_L2 + V(:,ii) .* W(:,ii);
    Ct_H1 = Ct_H1 + Vu(:,ii) .* Wu(:,ii);
end
Ct_L2 = Ct_L2 .* Cspeed;
Ct_H1 = Ct_H1 .* CspeedInv;

% H2 Energy terms
CuCuu = Cu(:,1) .* Cuu(:,1);
VuWu = Vu(:,1) .* Wu(:,1);   
VuWuu = Vu(:,1) .* Wuu(:,1);   
VuuWu = Vuu(:,1) .* Wu(:,1);   
VuuWuu = Vuu(:,1) .* Wuu(:,1);

for ii = 2:dSpace
    CuCuu = CuCuu + Cu(:,ii) .* Cuu(:,ii);
    VuWu = VuWu + Vu(:,ii) .* Wu(:,ii);
    VuWuu = VuWuu + Vu(:,ii) .* Wuu(:,ii);
    VuuWu = VuuWu + Vuu(:,ii) .* Wu(:,ii);
    VuuWuu = VuuWuu + Vuu(:,ii) .* Wuu(:,ii);
end

% Ct_H2 = CutCut .* CuCuu.^2 ./ Cspeed.^7 ...
%     - 2 * CutCuut .* CuCuu ./ Cspeed .^ 5 ...
%     + CuutCuut ./ Cspeed .^ 3;

CspeedInv2 = CspeedInv .* CspeedInv;
CspeedInv3 = CspeedInv2 .* CspeedInv;
CspeedInv5 = CspeedInv3 .* CspeedInv2;
CspeedInv7 = CspeedInv5 .* CspeedInv2;
Ct_H2 = VuWu .* CuCuu.^2 .* CspeedInv7 ...
    - (VuWuu + VuuWu) .* CuCuu .* CspeedInv5 ...
    + VuuWuu .* CspeedInv3;

a = splineData.a;
energyIntegrand = a(1) * Ct_L2 + a(2) * Ct_H1 + a(3) * Ct_H2;

%Compute final energy
G = sum(quadData.quadWeightsS .* energyIntegrand);

end

