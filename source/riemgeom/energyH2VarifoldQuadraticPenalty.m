%% energyH2Varifold
%
% Computes the sum
%   E(dPath) + lambda * dist(dPath(1), gamma.dEnd)^2
% with E(dPath) the Riemannian path energy and dist(.,.) the varifold
% distance between the endpoint of the path and the target curve
% transformed by gamma.
%
% Metric: G_c(h,h) = int a0<h,h> + (a1+a1n)<D_s h,D_s h> + ...
%                       (a1v-a1n)<D_s h,v>^2 + a2<D^2_s h,D^2_s h> ds
%
% If splineData.scaleInv = 1
%
% Metric: G_c(h,h) = int a0/l^3<h,h> + (a1+a1n)/l<D_s h,D_s h> + ...
%                       (a1v-a1n)/l<D_s h,v>^2 + a2*l<D^2_s h,D^2_s h> ds
%
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
function [E, dE] = energyH2VarifoldQuadraticPenalty(coeffs, params, splineData)

%% Extract parameters from splineData
dSpace = splineData.dSpace;
N = splineData.N;
Nt = splineData.Nt;
quadDataTensor = splineData.quadDataTensor;

a = splineData.a; a0 = a(1); a1 = a(2); a2 = a(3); a1v = a(4); a1n = a(5);

scaleInv = 0;
if ~isempty(splineData.scaleInv)
    scaleInv = splineData.scaleInv;
end

noQuadPointsS = splineData.quadData.noQuadPointsS;
noQuadPointsT = splineData.quadData.noQuadPointsT;
quadWeightsS = splineData.quadData.quadWeightsS;

%% Decode coeffs and params
d0 = params.d0;
dEnd = params.dEnd;
optTra = params.optTra;
optRot = params.optRot;
optScal = params.optScal;
lambda = params.lambda; %Weight of Varifold term in the energy functional

dPath = [ d0;
          reshape(coeffs(1:N*(Nt-1)*dSpace), [N*(Nt-1), dSpace]) ];
rho = coeffs(end-dSpace-1);
beta = coeffs(end-dSpace);
v = coeffs(end-dSpace+1:end);

%% Check whether path changes turning number
% This indicates that hanso wants to take a too large step length
checkTurningNumber = params.checkTurningNumber;
if checkTurningNumber
    change = changeTurningNumber(dPath, splineData);
    if change
        E =  inf;
        dE = inf;
        return;
    end
end

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
CutCu = sum(Cut.*Cu,2);
CutCu2 = CutCu.^2;

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
Ct_H1v = CutCu2.*CspeedInv3;

% H2 Energy terms
Ct_H2 = CutCut .* CuCuu2 .* CspeedInv7 ...
    - 2 * CutCuut .* CuCuu .* CspeedInv5 ...
    + CuutCuut .* CspeedInv3;

if scaleInv == 0
    Ct_L2H1H2 = a0*Ct_L2 +(a1+a1n)*Ct_H1 + (a1v-a1n)*Ct_H1v + a2*Ct_H2;
    energyIntegrand = Ct_L2H1H2;
else
    % Compute length-weighted metric
    ellVec = quadWeightsS'*reshape(Cspeed,noQuadPointsS,noQuadPointsT);
    ellInvVec = ellVec.^(-1);
    ellInv3Vec = ellVec.^(-3);
    
    ell = reshape( repmat( ellVec , noQuadPointsS, 1), [],1) ;
    ellInv = reshape( repmat( ellInvVec , noQuadPointsS, 1), [],1);
    ellInv3 = reshape( repmat( ellInv3Vec , noQuadPointsS, 1), [],1);
    
    Ct_L2H1H2 = a0*ellInv3.*Ct_L2 +(a1+a1n)*ellInv.*Ct_H1 ...
        + (a1v-a1n)*ellInv.*Ct_H1v + a2*ell.*Ct_H2;
    energyIntegrand = Ct_L2H1H2;
end
% Compute final energy
E = quadDataTensor.quadWeights' * energyIntegrand;


%% Compute Varifold energy
rotation = [ cos(beta), sin(beta); ...
             -sin(beta), cos(beta) ];

d2 = dEnd + ones([N, 1]) * v(:)';         
d2 = d2 * rotation;
d2 = rho*d2;


d1 = dPath(end-N+1:end,:);
[distVar, distGradd1, distGradd2] = varifoldDistanceSquared(d1, d2, splineData);

%% Compute the final energy
E = E + lambda*distVar;

%% Compute Gradient
if nargout > 1
    
    if scaleInv == 0
        term1 = a0*CspeedInv.*CtCt ...
            - (a1+a1n)*CspeedInv3.*CutCut ...
            - a2*7*CspeedInv9.*(CuCuu2).*CutCut ...
            + a2*10*CspeedInv7.*CuCuu.*CutCuut ...
            - a2*3*CspeedInv5.*CuutCuut ...
            - (a1v-a1n)*3*CspeedInv5.*CutCu2;
        term2 = a2*2*CspeedInv7.*CuCuu.*CutCut ...
            -a2*2*CspeedInv5.*CutCuut;
        term3 = 2*a0*Cspeed;
        term4 = 2*(a1+a1n)*CspeedInv + 2*a2*CspeedInv7.*CuCuu.^2;
        term5 = -a2*2*CspeedInv5.*CuCuu;
        term6 = 2*a2*CspeedInv3;
        term7 = 2*(a1v-a1n)*CspeedInv3.*CutCu;
        
        term1 = term1.*quadDataTensor.quadWeights;
        term2 = term2.*quadDataTensor.quadWeights;
        term3 = term3.*quadDataTensor.quadWeights;
        term4 = term4.*quadDataTensor.quadWeights;
        term5 = term5.*quadDataTensor.quadWeights;
        term6 = term6.*quadDataTensor.quadWeights;
        term7 = term7.*quadDataTensor.quadWeights;
        
        dE = zeros( splineData.N*splineData.Nt, splineData.dSpace);
        for kk = splineData.dSpace:-1:1
            termBu = (term1.*Cu(:,kk) + term2.*Cuu(:,kk) + term7.*Cut(:,kk))'*quadDataTensor.Bu;
            termBuu = (term2.*Cu(:,kk))'*quadDataTensor.Buu;
            termBt = (term3.*Ct(:,kk))'*quadDataTensor.Bt;
            termBut = (term4.*Cut(:,kk)+term5.*Cuut(:,kk)+term7.*Cu(:,kk))'*quadDataTensor.But;%last term here is new
            termBuut = (term6.*Cuut(:,kk)+term5.*Cut(:,kk))'*quadDataTensor.Buut;
            
            dE(:,kk) = termBu + termBuu + termBt + termBut + termBuut;
        end
  
    elseif scaleInv == 1
        ellInv2Vec = ellVec.^(-2);
        ellInv2 = reshape( repmat( ellInv2Vec , noQuadPointsS, 1), [],1);
        ellInv4Vec = ellVec.^(-4);
        ellInv4 = reshape( repmat( ellInv4Vec , noQuadPointsS, 1), [],1);
        
        Ct_L2_IntTheta = quadWeightsS'*reshape(Ct_L2,noQuadPointsS,noQuadPointsT);
        Ct_H1_IntTheta = quadWeightsS'*reshape(Ct_H1,noQuadPointsS,noQuadPointsT);
        Ct_H1v_IntTheta = quadWeightsS'*reshape(Ct_H1v,noQuadPointsS,noQuadPointsT);
        Ct_H2_IntTheta = quadWeightsS'*reshape(Ct_H2,noQuadPointsS,noQuadPointsT);
        
        Ct_L2_IntTheta = reshape( repmat( Ct_L2_IntTheta , noQuadPointsS, 1), [],1) ;
        Ct_H1_IntTheta = reshape( repmat( Ct_H1_IntTheta , noQuadPointsS, 1), [],1) ;
        Ct_H1v_IntTheta = reshape( repmat( Ct_H1v_IntTheta , noQuadPointsS, 1), [],1) ;
        Ct_H2_IntTheta = reshape( repmat( Ct_H2_IntTheta , noQuadPointsS, 1), [],1) ;
        
        term1 = a0*ellInv3.*CspeedInv.*CtCt ...
            - (a1+a1n)*ellInv.*CspeedInv3.*CutCut ...
            - a2*7*ell.*CspeedInv9.*(CuCuu2).*CutCut ...
            + a2*10*ell.*CspeedInv7.*CuCuu.*CutCuut ...
            - a2*3*ell.*CspeedInv5.*CuutCuut ...
            - (a1v-a1n)*3*ellInv.*CspeedInv5.*CutCu2 ...
            - a0*3*ellInv4.*Ct_L2_IntTheta.*CspeedInv ... %L-W L2 term
            - (a1+a1n)*ellInv2.*Ct_H1_IntTheta.*CspeedInv ...
            - (a1v-a1n)*ellInv2.*Ct_H1v_IntTheta.*CspeedInv ...
            + a2*Ct_H2_IntTheta.*CspeedInv;
        term2 = a2*2*ell.*CspeedInv7.*CuCuu.*CutCut ...
            -a2*2*ell.*CspeedInv5.*CutCuut;
        term3 = 2*a0*ellInv3.*Cspeed;
        term4 = 2*(a1+a1n)*ellInv.*CspeedInv + 2*a2*ell.*CspeedInv7.*CuCuu.^2;
        term5 = -a2*2*ell.*CspeedInv5.*CuCuu;
        term6 = 2*a2*ell.*CspeedInv3;
        term7 = 2*(a1v-a1n)*ellInv.*CspeedInv3.*CutCu;
        
        term1 = term1.*quadDataTensor.quadWeights;
        term2 = term2.*quadDataTensor.quadWeights;
        term3 = term3.*quadDataTensor.quadWeights;
        term4 = term4.*quadDataTensor.quadWeights;
        term5 = term5.*quadDataTensor.quadWeights;
        term6 = term6.*quadDataTensor.quadWeights;
        term7 = term7.*quadDataTensor.quadWeights;
        
        dE = zeros( splineData.N*splineData.Nt, splineData.dSpace);
        for kk = splineData.dSpace:-1:1
            termBu = (term1.*Cu(:,kk) + term2.*Cuu(:,kk) + term7.*Cut(:,kk))'*quadDataTensor.Bu;
            termBuu = (term2.*Cu(:,kk))'*quadDataTensor.Buu;
            termBt = (term3.*Ct(:,kk))'*quadDataTensor.Bt;
            termBut = (term4.*Cut(:,kk)+term5.*Cuut(:,kk)+term7.*Cu(:,kk))'*quadDataTensor.But;
            termBuut = (term6.*Cuut(:,kk)+term5.*Cut(:,kk))'*quadDataTensor.Buut;
            
            dE(:,kk) = termBu + termBuu + termBt + termBut + termBuut;
        end
    end
    dEdc1 = dE(end-splineData.N+1:end,:);

    dE = dE(splineData.N+1:end-splineData.N,:);
    
    % Compute Varifold Gradient
    % distGradd1, distGradd2 has been computed above
    % [~, distGradd1] = varifoldDistanceSquared(d1, d2, splineData);
    %if optRot || optTra || optScal
    %    [~, distGradd2] = varifoldDistanceSquared(d2, d1, splineData);
    %end
    
    % Compute gradient w.r.t beta and v
    % Observe that F(rho,beta,v) = || d0 - rho R_beta (d1 + v) ||^2 

    if optRot
        rotationDer = [ -sin(beta), cos(beta); -cos(beta), -sin(beta) ]; 
        distVarGradBeta = rho*sum(sum(distGradd2 .*((dEnd+ones([N, 1]) * v(:)')*rotationDer)));
    else
        distVarGradBeta = 0; 
    end
    
    if optScal 
        distVarGradRho = sum(sum(distGradd2 .*((dEnd+ones([N, 1]) * v(:)')*rotation)));
    else
        distVarGradRho = 0; 
    end
    
    if optTra
        distVarGradV = rho*rotation*sum(distGradd2,1 )';%rho.* 
    else
        distVarGradV = [0; 0];
    end
    
    % Update the endpoint part with the derivatives of varifold distance
    dE = [ dE; dEdc1 + lambda*distGradd1];
  
    % Attach gradients wrt rho, beta and v
    dE = [ dE(:);lambda*distVarGradRho; lambda*distVarGradBeta; lambda*distVarGradV ]; 
end

end

