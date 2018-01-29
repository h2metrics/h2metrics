%% pathRiemH2Energy
%
% Computes the Riemannian energy of a spline path using the formula
%   \int_0^{2pi} a(1) * <h, h> + (a(2)+a(5)) * <D_s h, D_s h> + ...
%                (a(4)-a(5)) * <D_s h, v>^2 + a(3) * <D^2_s h, D^2_s h> ds
%
% Input
%   dPath
%       The spline path
%   splineData
%       General information about the splines used.
%
% Optional inputs
%   a
%       Coefficients for the inner product
%
% Output
%   G
%       The energy of the path
%   comp
%       The three components of the norm separately
%           [ L2, H1, H2, H1v, H1n ]
%       We have the identity
%           G = L2 + H1 + H2 + H1v + H1n
%
% Notes
%   The order of precedence for the constants are as follows
%     -) Optional parameter 'a'
%     -) splineData.a
%     -) a = [1 0 1 0 0]
%
function [E, comp] = pathRiemH2Energy( dPath, splineData, varargin )
                                   
% Set constants to be used in metric
a = [1 0 1 0 0];
if ~isempty(splineData.a)
    a = splineData.a;
end

% Parse optional inputs
p = inputParser;
addParameter(p, 'a', a);
parse(p, varargin{:});
a = p.Results.a;

quadData = splineData.quadData;
quadDataTensor = splineData.quadDataTensor;

%% Evaluate path at quadrature sites
Cu = quadDataTensor.Bu * dPath;
Ct = quadDataTensor.Bt * dPath;
Cut = quadDataTensor.But * dPath;
Cuu = quadDataTensor.Buu * dPath;
Cuut = quadDataTensor.Buut * dPath;

%% Calculate terms of the energy
CuCuu = sum(Cu.*Cuu,2);
CutCut = sum(Cut.*Cut,2);
CutCuut = sum(Cut.*Cuut,2);
CuutCuut = sum(Cuut.*Cuut,2);
CutCu = sum(Cut.*Cu,2);
CutCu2 = CutCu.^2;

Cspeed = sum( Cu.^2 , 2).^(1/2);

CspeedInv = 1./Cspeed;
CspeedInv2 = CspeedInv .* CspeedInv;
CspeedInv3 = CspeedInv2 .* CspeedInv;
CspeedInv5 = CspeedInv3 .* CspeedInv2;
CspeedInv7 = CspeedInv5 .* CspeedInv2;

Ct_L2 = sum( Ct .* Ct, 2) .* Cspeed;
Ct_H1 = sum( Cut .* Cut, 2) .* CspeedInv;
Ct_H1v = CutCu2 .* CspeedInv3;
Ct_H2 = CutCut .* CuCuu.^2 .* CspeedInv7 ...
    - 2 * CutCuut .* CuCuu .* CspeedInv5 ...
    + CuutCuut .* CspeedInv3;

%% Now integrate
scaleInv = splineData.scaleInv;
if ~scaleInv
    L2 = sum(Ct_L2 .* quadDataTensor.quadWeights);
    H1 = sum(Ct_H1 .* quadDataTensor.quadWeights);
    H2 = sum(Ct_H2 .* quadDataTensor.quadWeights);
    H1v = sum(Ct_H1v .* quadDataTensor.quadWeights);
    H1n = H1 - H1v;
else
    noQuadPointsS = splineData.quadData.noQuadPointsS;
    noQuadPointsT = splineData.quadData.noQuadPointsT;
    quadWeightsS = quadData.quadWeightsS;
    
    ellVec = quadWeightsS' * reshape(Cspeed, noQuadPointsS, noQuadPointsT);
    ellInvVec = ellVec.^(-1);
    ellInv3Vec = ellVec.^(-3);
    
    ell = reshape( repmat( ellVec , noQuadPointsS, 1), [],1) ;
    ellInv = reshape( repmat( ellInvVec , noQuadPointsS, 1), [],1);
    ellInv3 = reshape( repmat( ellInv3Vec , noQuadPointsS, 1), [],1);
    
    L2 = sum(ellInv3 .* Ct_L2 .* quadDataTensor.quadWeights);
    H1 = sum(ellInv .* Ct_H1 .* quadDataTensor.quadWeights);
    H2 = sum(ell .* Ct_H2 .* quadDataTensor.quadWeights);
    H1v = sum(ellInv .* Ct_H1v .* quadDataTensor.quadWeights);
    H1n = H1 - H1v;
end

E = a(1) * L2 + a(2) * H1 + a(3) * H2 + a(4) * H1v + a(5) * H1n;
comp = [a(1)*L2, a(2)*H1, a(3)*H2, a(4)*H1v, a(5)*H1n];

end