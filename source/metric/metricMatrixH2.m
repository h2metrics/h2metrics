%% metricMatrixH2
%
% Function computes the metric matrix at the curve d. The metric matrix
% corresponds to the flattening of d as
%   [d(:,1), d(:,2), ..., d(:,dSpace)]
%
% Input
%   d
%       The curve. Has dimensions [N, dSpace]
%   splineData
%       splineData describing the curve
%
% Output
%   G
%       Matrix of the Riemannian metric. Has dimensions 
%         [N*dSpace, N*dSpace]
%
function [ G ] = metricMatrixH2( d, splineData )

a = splineData.a;

N = splineData.N;
nS = splineData.nS;
dSpace = splineData.dSpace;
quadData = splineData.quadData;
noQuadSites = quadData.noQuadPointsS;
quadWeights = quadData.quadWeightsS;

scaleInv = 0;
if ~isempty(splineData.scaleInv)
    scaleInv = splineData.scaleInv;
end

B = quadData.B_S;
Bu = quadData.Bu_S;
Buu = quadData.Buu_S;

Cu = quadData.Bu_S * d;
Cuu = quadData.Buu_S * d;

CuCuu = sum(Cu.*Cuu,2);
CuCuu2 = CuCuu .* CuCuu;

Cspeed = sum( Cu.^2 , 2).^(1/2);
CspeedInv = 1 ./ Cspeed;
CspeedInv2 = CspeedInv .* CspeedInv;
CspeedInv3 = CspeedInv .* CspeedInv2;
CspeedInv5 = CspeedInv3 .* CspeedInv2;
CspeedInv7 = CspeedInv5 .* CspeedInv2;

% Collect the terms in the metric
G = spalloc(N*dSpace, N*dSpace, N*(nS+1)*dSpace);

weight_L2 = Cspeed .* quadWeights;
sparseW_L2 = sparse(1:noQuadSites, 1:noQuadSites, weight_L2);

weight_H1 = CspeedInv .* quadWeights;
sparseW_H1 = sparse(1:noQuadSites, 1:noQuadSites, weight_H1);

% <D_s h,v><D_s k, v>|c'| = 1/|c'|^3 <h',c'> <k',c'>
sparseW_H1v = sparse(1:noQuadSites, 1:noQuadSites, weight_H1);
sparseCu = sparse(1:noQuadSites*dSpace,repmat(1:noQuadSites,1,dSpace), Cu(:) );
sparseCuBu = spalloc(  noQuadSites, dSpace*N, dSpace*nnz(Bu) );

sparseCuBu = cell(2,1);
for ii = 1:dSpace
    sparseCuii = sparse(1:noQuadSites,1:noQuadSites, Cu(:,ii) );
    sparseCuBu{ii} = sparseCuii*Bu;    
%     sparseCuBu( 1:noQuadSites,(1+(ii-1)*N):ii*N )  = sparseCuii*Bu;    
end

% sparse( 1:2*noQuadSites , 

% reshape( sparseCuBu', 2

weight_H2_1 = CspeedInv7 .* CuCuu2 .* quadWeights;
weight_H2_2 = -CspeedInv5 .* CuCuu .* quadWeights;
weight_H2_3 = CspeedInv3 .* quadWeights;
sparseW_H2_1 = sparse(1:noQuadSites, 1:noQuadSites, weight_H2_1);
sparseW_H2_2 = sparse(1:noQuadSites, 1:noQuadSites, weight_H2_2);
sparseW_H2_3 = sparse(1:noQuadSites, 1:noQuadSites, weight_H2_3);

if scaleInv
    ell = quadWeights' * Cspeed; %Length of c
    %Update coeffecients with the length weights
    a(1) = a(1)/ell^3;
    a(2) = a(2)/ell;
    a(3) = a(3)*ell;
    a(4) = a(4)/ell;
    a(5) = a(5)/ell;
end

G_L2 = a(1) * B'*sparseW_L2*B;
G_H1 = a(2) * Bu'*sparseW_H1*Bu;
G_H2 = a(3) * (Bu'*sparseW_H2_1*Bu + Bu'*sparseW_H2_2*Buu + ...
               Buu'*sparseW_H2_2*Bu + Buu'*sparseW_H2_3*Buu);        

% Check the space dimension, if it is not 2, then the G_H1v and G_H1n terms
% have to be updated.
if dSpace > 2
    error('This function does not support 3D yet, update normal projection')
end
sparseTanBu = [sparseCuBu{1},sparseCuBu{2}];
sparseNormalBu = [-sparseCuBu{2},sparseCuBu{1}];
G_H1v = a(4)*sparseTanBu'*sparseW_H1*sparseTanBu; %Full size of G           
G_H1n = a(5)*sparseNormalBu'*sparseW_H1*sparseNormalBu; %Full size of G           

% Add L2,H1 and H2 terms to metric on the diagonal
G = zeros(N*dSpace, N*dSpace);
for kk = 1:dSpace 
    G( (kk-1)*N+1:kk*N, (kk-1)*N+1:kk*N ) = G_L2 + G_H1 + G_H2;
end

% Add the (a,b) terms
G = G + G_H1v + G_H1n;

end


