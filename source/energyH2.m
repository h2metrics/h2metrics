%% energyH2Diff
%
% Computes the H2 energy of a path; uses preconditioner. Function is used
% by geodesicBvpNoGroups
%
% Input
%   coeffs
%       Path of curves; matrix of dimension [N*(Nt-2), dSpace]
%   params
%       Start and end point of curve, preconditioner
%
% Output
%   E
%       Energy of the path.
%   dE
%       Derivative of the energy
%   H
%       Hessian of the energy
%
function [E, dE, H] = energyH2( coeffs, params, splineData)

%% Extract parameters
dSpace = splineData.dSpace;
N = splineData.N;
Nt = splineData.Nt;
quadDataTensor = splineData.quadDataTensor;

a = splineData.a; a0 = a(1); a1 = a(2); a2 = a(3);

noControlPoints = splineData.N*splineData.Nt;
noControls = noControlPoints*splineData.dSpace;
noVariables = splineData.N*(splineData.Nt-2); % BVP
noQuadSites = length( quadDataTensor.quadWeights );

%% Decode coeffs and params
d0 = params.d0;
d1 = params.d1;

PinvBlock = params.PinvBlock; % Preconditioning
coeffs = PinvBlock * coeffs;

dPath = [ d0; 
          reshape(coeffs(1:N*(Nt-2)*dSpace), [N*(Nt-2), dSpace]); 
          d1 ];

dPath2 = dPath;

%% Evaluate path at quadrature sites
Cu = quadDataTensor.Bu*dPath2;
Ct = quadDataTensor.Bt*dPath2;
Cut = quadDataTensor.But*dPath2;
Cuu = quadDataTensor.Buu*dPath2;
Cuut = quadDataTensor.Buut*dPath2;

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
    
    dE = dE(splineData.N+1:end-splineData.N,:);
    
    dE = dE(:);
    
    %% Add preconditioner
    dE = PinvBlock' * dE;
    
end

%% Compute Hessian
if nargout > 2
    CspeedInv11 = CspeedInv9.*CspeedInv2;
    
    weight1 = -a0*CspeedInv3.*CtCt + 3*a1*CspeedInv5.*CutCut ...
        + 63*a2*CspeedInv11.*CuCuu2.*CutCut ...
        - 70*a2*CspeedInv9.*CuCuu.*CutCuut ...
        + 15*a2*CspeedInv7.*CuutCuut;
    weight2 = -14*a2*CspeedInv9.*CuCuu.*CutCut + 10*a2*CspeedInv7.*CutCuut;
    weight3 = 2*a2*CspeedInv7.*CutCut;
    weight4 = 2*a0*CspeedInv;
    weight5 = -2*a1*CspeedInv3-14*a2*CspeedInv9.*CuCuu2;
    weight6 = 10*a2*CspeedInv7.*CuCuu;   
    weight7 = 4*a2*CspeedInv7.*CuCuu; 
    weight8 = -2*a2*CspeedInv5;
    weight9 = -6*a2*CspeedInv5;
   
    weight1 = weight1.*quadDataTensor.quadWeights;
    weight2 = weight2.*quadDataTensor.quadWeights;
    weight3 = weight3.*quadDataTensor.quadWeights;
    weight4 = weight4.*quadDataTensor.quadWeights;
    weight5 = weight5.*quadDataTensor.quadWeights;
    weight6 = weight6.*quadDataTensor.quadWeights;
    weight7 = weight7.*quadDataTensor.quadWeights;
    weight8 = weight8.*quadDataTensor.quadWeights;
    weight9 = weight9.*quadDataTensor.quadWeights;
    
    %Allocate sparse Hessian
    nnzHess = quadDataTensor.nnz;% nnz( quadDataTensor.B'*quadDataTensor.B ); %TODO: Cheaper 
    Hess = spalloc( noControls,noControls, nnzHess*splineData.dSpace^2);
    Hess2 = spalloc( noControls,noControls, nnzHess*splineData.dSpace);
    
    for kk_h = 1:splineData.dSpace
       for kk_u = 1:splineData.dSpace
           %Symmetric
           %Hu'*(...)*Uu
           weightHuUu = weight1.*Cu(:,kk_h).*Cu(:,kk_u) + ...
               weight2.*Cuu(:,kk_h).*Cu(:,kk_u) + ...
               weight2.*Cu(:,kk_h).*Cuu(:,kk_u) + ...
               weight3.*Cuu(:,kk_h).*Cuu(:,kk_u);
           sparseWHuUu = sparse(1:noQuadSites,1:noQuadSites,...
               weightHuUu );

%            termHuUu = quadDataTensor.BuTr*sparseWHuUu*quadDataTensor.Bu;
         
           %Huu'*(...)*Uuu
           weightHuuUuu = weight3.*Cu(:,kk_h).*Cu(:,kk_u);
           sparseWHuuUuu = sparse(1:noQuadSites,1:noQuadSites,...
               weightHuuUuu );
%            termHuuUuu = quadDataTensor.BuuTr*sparseWHuuUuu*quadDataTensor.Buu;
         
           %Mixed
           %Huu'*(...)*Uu
           weightHuuUu = weight2.*Cu(:,kk_h).*Cu(:,kk_u) + ...
               weight3.*Cu(:,kk_h).*Cuu(:,kk_u);
           sparseWHuuUu = sparse(1:noQuadSites,1:noQuadSites,...
               weightHuuUu );
%            termHuuUu = quadDataTensor.BuuTr*sparseWHuuUu*quadDataTensor.Bu;
           
           %Ht'*(...)*Uu
           weightHtUu = weight4.*Ct(:,kk_h).*Cu(:,kk_u);
           sparseWHtUu = sparse(1:noQuadSites,1:noQuadSites,...
               weightHtUu );
%            termHtUu = quadDataTensor.BtTr*sparseWHtUu*quadDataTensor.Bu;
           
           %Hut'*(...)*Uu
           weightHutUu = weight5.*Cut(:,kk_h).*Cu(:,kk_u) + ...
               weight6.*Cuut(:,kk_h).*Cu(:,kk_u) + ...
               weight7.*Cut(:,kk_h).*Cuu(:,kk_u) + ...
               weight8.*Cuut(:,kk_h).*Cuu(:,kk_u);
           sparseWHutUu = sparse(1:noQuadSites,1:noQuadSites,...
               weightHutUu );
%            termHutUu = quadDataTensor.ButTr*sparseWHutUu*quadDataTensor.Bu;
           
           %Huut'*(...)*Uu
           weightHuutUu = weight6.*Cut(:,kk_h).*Cu(:,kk_u) + ...
               weight9.*Cuut(:,kk_h).*Cu(:,kk_u) + ...
               weight8.*Cut(:,kk_h).*Cuu(:,kk_u);
           sparseWHuutUu = sparse(1:noQuadSites,1:noQuadSites,...
               weightHuutUu );
%            termHuutUu = quadDataTensor.BuutTr*sparseWHuutUu*quadDataTensor.Bu;
           
           %Hut'*(...)*Uuu
           weightHutUuu = weight7.*Cut(:,kk_h).*Cu(:,kk_u) + ...
               weight8.*Cuut(:,kk_h).*Cu(:,kk_u);
           sparseWHutUuu = sparse(1:noQuadSites,1:noQuadSites,...
               weightHutUuu );
%            termHutUuu = quadDataTensor.ButTr*sparseWHutUuu*quadDataTensor.Buu;
           
           %Huut'*(...)*Uuu
           weightHuutUuu = weight8.*Cut(:,kk_h).*Cu(:,kk_u);
           sparseWHuutUuu = sparse(1:noQuadSites,1:noQuadSites,...
               weightHuutUuu );
%            termHuutUuu = quadDataTensor.BuutTr*sparseWHuutUuu*quadDataTensor.Buu;
               
           totalterms = (1/2*quadDataTensor.BuTr*sparseWHuUu + ...
               quadDataTensor.BuuTr*sparseWHuuUu + ...
               quadDataTensor.BtTr*sparseWHtUu + ...
               quadDataTensor.ButTr*sparseWHutUu + ...
               quadDataTensor.BuutTr*sparseWHuutUu)*quadDataTensor.Bu + ...
               (1/2*quadDataTensor.BuuTr*sparseWHuuUuu + ...
               quadDataTensor.ButTr*sparseWHutUuu + ...
               quadDataTensor.BuutTr*sparseWHuutUuu)*quadDataTensor.Buu;
           
           Hess( (kk_h-1)*noControlPoints+1:kk_h*noControlPoints,...
               (kk_u-1)*noControlPoints+1:kk_u*noControlPoints) = totalterms;
%            temp = termHuuUu + termHtUu + termHutUu + termHuutUu + ...
%                termHutUuu + termHuutUuu;
%            Hess( (kk_h-1)*noControlPoints+1:kk_h*noControlPoints,...
%                (kk_u-1)*noControlPoints+1:kk_u*noControlPoints) ...
%                = 1/2*(termHuUu + termHuuUuu) + temp;
       end
    end
    Hess = Hess + Hess';
    
    %Symmetric
    sparseWHuUu = sparse(1:noQuadSites,1:noQuadSites, term1 );
    termHuUu = quadDataTensor.BuTr*sparseWHuUu*quadDataTensor.Bu;
        
    sparseWHtUt = sparse(1:noQuadSites,1:noQuadSites, term3 );
    termHtUt = quadDataTensor.BtTr*sparseWHtUt*quadDataTensor.Bt;
    
    sparseWHutUut = sparse(1:noQuadSites,1:noQuadSites, term4 );
    termHutUut = quadDataTensor.ButTr*sparseWHutUut*quadDataTensor.But;
    
    sparseWHuutUuut = sparse(1:noQuadSites,1:noQuadSites, term6 );
    termHuutUuut = quadDataTensor.BuutTr*sparseWHuutUuut*quadDataTensor.Buut;
    
    %Mixed
    sparseWHuuUu = sparse(1:noQuadSites,1:noQuadSites, term2 );
    termHuuUu = quadDataTensor.BuuTr*sparseWHuuUu*quadDataTensor.Bu;
    
    sparseWHuutUut = sparse(1:noQuadSites,1:noQuadSites, term5 );
    termHuutUut = quadDataTensor.BuutTr*sparseWHuutUut*quadDataTensor.But;
    
    temp = termHuuUu + termHuutUut;
    temp2 = termHuUu + termHtUt + termHutUut + termHuutUuut + temp + temp';
    for kk = 1:splineData.dSpace
        Hess2( (kk-1)*noControlPoints+1:kk*noControlPoints,...
               (kk-1)*noControlPoints+1:kk*noControlPoints) ...
               = temp2;
    end
    Hess = Hess + Hess2;
    
    %Extract indices for interior controls in H, upper left block
    intIndices  = zeros( 1,noVariables );
    for kk = 1:splineData.dSpace
        intIndices( (kk-1)*noVariables+1:kk*noVariables) = ...
            (kk-1)*noControlPoints+[splineData.N+1:noControlPoints-splineData.N];
    end
        
    HessInt = Hess( intIndices, intIndices);
    
%     for ii = dSpace:-1:1
%         for jj = Nt-2:-1:1
%             for kk = dSpace:-1:1
%                 for ll = Nt-2:-1:1
%                     HessInt((ii-1)*N*(Nt-2)+(jj-1)*N+1:...
%                             (ii-1)*N*(Nt-2)+jj*N, ...
%                             (kk-1)*N*(Nt-2)+(ll-1)*N+1:...
%                             (kk-1)*N*(Nt-2)+ll*N) = ...
%                         Pinv' * HessInt((ii-1)*N*(Nt-2)+(jj-1)*N+1:...
%                                         (ii-1)*N*(Nt-2)+jj*N, ...
%                                         (kk-1)*N*(Nt-2)+(ll-1)*N+1:...
%                                         (kk-1)*N*(Nt-2)+ll*N) * Pinv;
%                 end
%             end
%         end
%     end
    
    H = HessInt;   
    
    % NOTE: Assuming symmetry of PinvBlock
    H = PinvBlock * H * PinvBlock;
end



end

