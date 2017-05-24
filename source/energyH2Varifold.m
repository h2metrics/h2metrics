%% energyH2Varifold
%
% Computes the sum
%   E(dPath) + lambda * dist(dPath(1), gamma.dEnd)^2
% with E(dPath) the Riemannian path energy and dist(.,.) the varifold
% distance between the endpoint of the path and the target curve
% transformed by gamma.
%
% Input
%   coeffs
%       Path of curves, rotation and translation vector
%       Vector of length N*(Nt-1)*dSpace+dSpace+1
%   params
%       Start point of path, target curve, toggles for translation and
%       rotation
%   splineData
%       General information about the splines used.
%
% Output
%   E
%       Energy of the path.
%   dE
%       Derivative of the energy
%
function [E, dE] = energyH2Varifold(coeffs, params, splineData)

%% Extract parameters from splineData
dSpace = splineData.dSpace;
N = splineData.N;
Nt = splineData.Nt;

lambda = splineData.lambda; %Weight of Varifold term in the energy functional
quadDataTensor = splineData.quadDataTensor;
a = splineData.a; a0 = a(1); a1 = a(2); a2 = a(3);

%% Decode coeffs and params
d0 = params.d0;
dEnd = params.dEnd;
optTra = params.optTra;
optRot = params.optRot;

dPath = [ d0;
          reshape(coeffs(1:N*(Nt-1)*dSpace), [N*(Nt-1), dSpace]) ];
beta = coeffs(end-dSpace);
v = coeffs(end-dSpace+1:end);

%% Evaluate path at quadrature sites
Cu = quadDataTensor.Bu*dPath;
Ct = quadDataTensor.Bt*dPath;
Cut = quadDataTensor.But*dPath;
Cuu = quadDataTensor.Buu*dPath;
Cuut = quadDataTensor.Buut*dPath;

CuCuu = sum(Cu.*Cuu,2);
CuCuu2 = CuCuu.^2;
CtCt = sum(Ct.*Ct,2);
CutCut = sum(Cut.*Cut,2);
CutCuut = sum(Cut.*Cuut,2);
CuutCuut = sum(Cuut.*Cuut,2);

Cspeed = sum( Cu.^2 , 2).^(1/2);

CspeedInv = 1./Cspeed;
CspeedInv2 = CspeedInv .* CspeedInv;
CspeedInv3 = CspeedInv2 .* CspeedInv;
CspeedInv5 = CspeedInv3 .* CspeedInv2;
CspeedInv7 = CspeedInv5 .* CspeedInv2;
CspeedInv9 = CspeedInv7 .* CspeedInv2;

%% Energy of the path
% L2 Energy terms
Ct_L2 = CtCt.*Cspeed;

% H1 Energy terms
Ct_H1 = CspeedInv.*CutCut;

% H2 Energy terms
Ct_H2 = CutCut .* CuCuu2 .* CspeedInv7 ...
    - 2 * CutCuut .* CuCuu .* CspeedInv5 ...
    + CuutCuut .* CspeedInv3;

Ct_L2H1H2 = a(1)*Ct_L2 +a(2)*Ct_H1+ a(3)*Ct_H2;
energyIntegrand = Ct_L2H1H2;

% Compute final energy
E = quadDataTensor.quadWeights' * energyIntegrand;

%% Compute Varifold energy
rotation = [ cos(beta), sin(beta); ...
             -sin(beta), cos(beta) ];
d2 = dEnd * rotation;

d2 = d2 + ones([N, 1]) * v(:)';

d1 = dPath(end-N+1:end,:);
distVar = varifoldDistanceSquared(d1, d2, splineData);

%% Compute the final energy
E = E + lambda*distVar;

%% Compute Gradient
if nargout > 1
    
    term1 = a0*CspeedInv.*CtCt ...
        - a1*CspeedInv3.*CutCut ...
        - a2*7*CspeedInv9.*(CuCuu2).*CutCut ...
        + a2*10*CspeedInv7.*CuCuu.*CutCuut ...
        - a2*3*CspeedInv5.*CuutCuut;
    term2 = a2*2*CspeedInv7.*CuCuu.*CutCut ...
        -a2*2*CspeedInv5.*CutCuut;
    term3 = 2*a0*Cspeed;
    term4 = 2*a1*CspeedInv + 2*a2*CspeedInv7.*CuCuu.^2;
    term5 = -a2*2*CspeedInv5.*CuCuu;
    term6 = 2*a2*CspeedInv3;
    
    term1 = term1.*quadDataTensor.quadWeights;
    term2 = term2.*quadDataTensor.quadWeights;
    term3 = term3.*quadDataTensor.quadWeights;
    term4 = term4.*quadDataTensor.quadWeights;
    term5 = term5.*quadDataTensor.quadWeights;
    term6 = term6.*quadDataTensor.quadWeights;
    
    dE = zeros( splineData.N*splineData.Nt, splineData.dSpace);
    for kk = splineData.dSpace:-1:1
        term1test = (term1.*Cu(:,kk) + term2.*Cuu(:,kk))'*quadDataTensor.Bu;
        term2test = (term2.*Cu(:,kk))'*quadDataTensor.Buu;
        term3test = (term3.*Ct(:,kk))'*quadDataTensor.Bt;
        term4test = (term4.*Cut(:,kk)+term5.*Cuut(:,kk))'*quadDataTensor.But;
        term5test = (term6.*Cuut(:,kk)+term5.*Cut(:,kk))'*quadDataTensor.Buut;

        dE(:,kk) = term1test + term2test + term3test + term4test + term5test;
    end
    
    dEdc1 = dE(end-splineData.N+1:end,:);

    dE = dE(splineData.N+1:end-splineData.N,:);
    
    % Compute Varifold Gradient
    [~, distGradd1] = varifoldDistanceSquared(d1, d2, splineData);
    if optRot || optTra
        [~, distGradd2] = varifoldDistanceSquared(d2, d1, splineData);
    end
    
    % Compute gradient w.r.t beta and v
    % Observe that F(beta,v) = || d0 - (R_beta(d1) + v) ||^2 
    %                        = || R_{-beta}(d0 - v) - d1 ||^2  
    if optRot
        rotationDer = [ -sin(beta), cos(beta); -cos(beta), -sin(beta) ]; 
        distVarGradBeta = sum(sum(distGradd2 .*(dEnd*rotationDer)));
    else
        distVarGradBeta = 0; 
    end
    if optTra
        distVarGradV = sum( distGradd2,1 )'; 
    else
        distVarGradV = [0; 0];
    end
    
    % Update the endpoint part with the derivatives of varifold distance
    dE = [ dE; dEdc1 + lambda*distGradd1];
  
    % Attach gradients wrt beta and v
    dE = [ dE(:); lambda*distVarGradBeta; lambda*distVarGradV ]; 
end

end

