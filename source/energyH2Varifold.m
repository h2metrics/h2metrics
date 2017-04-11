%% energyH2DiffVarifold
%
% Computes the energy of dPath after some reparametrizations.
%
% Input
%   dPath
%       Path of curves; matrix of dimension [N*Nt, dSpace]
%
%   d
%       Curve that end points is compared to by the varifold distance
%   splineData
%       General information about the splines used.
%
% Output
%   E1
%       Energy of the path.
%   dE
%       Derivative of the energy
%   H
%       Hessian of the energy
%
function [E,dE] = energyH2Varifold( dPath, d, splineData, varargin)

%% Optional parameters
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1; 
    end
end

%% Extract parameters
lambda = splineData.lambda; %Weight of Varifold term in the energy functional
dSpace = splineData.dSpace;
N = splineData.N;
quadData = splineData.quadData;
quadDataTensor = splineData.quadDataTensor;

a = splineData.a; a0 = a(1); a1 = a(2); a2 = a(3);
noInterpolPoints = size(quadData.B_interpolPhi,1);

noControlPoints = splineData.N*splineData.Nt;
noControls = noControlPoints*splineData.dSpace;
noVariables = splineData.N*(splineData.Nt-1); %BVP
noQuadSites = length( quadDataTensor.quadWeights );

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
% noEvalPoints = 100;
% 
% % End point curve
% d1 = dPath(end-N+1:end,:);
% 
% % Parameters for varifold code
% objfun.kernel_geom='gaussian';
% objfun.kernel_size_geom=splineData.varifoldKernelSize;
% objfun.kernel_grass='binet';
% 
% % Evaluate the curve
% evalS = linspace(0,2*pi,noEvalPoints+1)'; 
% evalS = evalS(1:end-1); % Remove the last repeated point
% evald1 = evalCurve(evalS,d1,splineData);
% evalC = evalCurve(evalS,d,splineData); %Matching curve
% 
% % Create connectivity matrix
% G = [(1:length( evalS ))', circshift( (1:length( evalS ))' , -1)];
% 
% d1vari = struct( 'x', evald1, 'G', G );
% Cvari = struct( 'x', evalC, 'G', G );
% 
% distVar = varifoldnorm(d1vari,Cvari,objfun);

d1 = dPath(end-N+1:end,:);
distVar = varifoldDistanceSquared(d1, d, splineData);

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
    
    Edc1 = dE(end-splineData.N+1:end,:);

    dE = dE(splineData.N+1:end-splineData.N,:);
    
    % Compute Varifold Gradient
%     distVarGradPts = dvarifoldnorm(d1vari,Cvari,objfun);
%     
%     % dist: d -Ev-> B*d -varNorm-> distVar
%     % grad(dist) = (varNorm'*B)
%     
%     distVarGrad = (distVarGradPts'*splineData.quadData.B_interpolS)';
    
    [~, distVarGrad] = varifoldDistanceSquared(d1, d, splineData);
    
    %Collect all dE terms
    dE = [dE; Edc1 + lambda*distVarGrad];
    
end

% %% Compute Hessian
% if nargout > 2
%     CspeedInv11 = CspeedInv9.*CspeedInv2;
%     
%     weight1 = -a0*CspeedInv3.*CtCt + 3*a1*CspeedInv5.*CutCut ...
%         + 63*a2*CspeedInv11.*CuCuu2.*CutCut ...
%         - 70*a2*CspeedInv9.*CuCuu.*CutCuut ...
%         + 15*a2*CspeedInv7.*CuutCuut;
%     weight2 = -14*a2*CspeedInv9.*CuCuu.*CutCut + 10*a2*CspeedInv7.*CutCuut;
%     weight3 = 2*a2*CspeedInv7.*CutCut;
%     weight4 = 2*a0*CspeedInv;
%     weight5 = -2*a1*CspeedInv3-14*a2*CspeedInv9.*CuCuu2;
%     weight6 = 10*a2*CspeedInv7.*CuCuu;   
%     weight7 = 4*a2*CspeedInv7.*CuCuu; 
%     weight8 = -2*a2*CspeedInv5;
%     weight9 = -6*a2*CspeedInv5;
%    
%     weight1 = weight1.*quadDataTensor.quadWeights;
%     weight2 = weight2.*quadDataTensor.quadWeights;
%     weight3 = weight3.*quadDataTensor.quadWeights;
%     weight4 = weight4.*quadDataTensor.quadWeights;
%     weight5 = weight5.*quadDataTensor.quadWeights;
%     weight6 = weight6.*quadDataTensor.quadWeights;
%     weight7 = weight7.*quadDataTensor.quadWeights;
%     weight8 = weight8.*quadDataTensor.quadWeights;
%     weight9 = weight9.*quadDataTensor.quadWeights;
%     
%     %Allocate sparse Hessian
%     nnzHess = quadDataTensor.nnz;% nnz( quadDataTensor.B'*quadDataTensor.B ); %TODO: Cheaper 
%     Hess = spalloc( noControls,noControls, nnzHess*splineData.dSpace^2);
%     Hess2 = spalloc( noControls,noControls, nnzHess*splineData.dSpace);
%     
%     for kk_h = 1:splineData.dSpace
%        for kk_u = 1:splineData.dSpace
%            %Symmetric
%            %Hu'*(...)*Uu
%            weightHuUu = weight1.*Cu(:,kk_h).*Cu(:,kk_u) + ...
%                weight2.*Cuu(:,kk_h).*Cu(:,kk_u) + ...
%                weight2.*Cu(:,kk_h).*Cuu(:,kk_u) + ...
%                weight3.*Cuu(:,kk_h).*Cuu(:,kk_u);
%            sparseWHuUu = sparse(1:noQuadSites,1:noQuadSites,...
%                weightHuUu );
% 
% %            termHuUu = quadDataTensor.BuTr*sparseWHuUu*quadDataTensor.Bu;
%          
%            %Huu'*(...)*Uuu
%            weightHuuUuu = weight3.*Cu(:,kk_h).*Cu(:,kk_u);
%            sparseWHuuUuu = sparse(1:noQuadSites,1:noQuadSites,...
%                weightHuuUuu );
% %            termHuuUuu = quadDataTensor.BuuTr*sparseWHuuUuu*quadDataTensor.Buu;
%          
%            %Mixed
%            %Huu'*(...)*Uu
%            weightHuuUu = weight2.*Cu(:,kk_h).*Cu(:,kk_u) + ...
%                weight3.*Cu(:,kk_h).*Cuu(:,kk_u);
%            sparseWHuuUu = sparse(1:noQuadSites,1:noQuadSites,...
%                weightHuuUu );
% %            termHuuUu = quadDataTensor.BuuTr*sparseWHuuUu*quadDataTensor.Bu;
%            
%            %Ht'*(...)*Uu
%            weightHtUu = weight4.*Ct(:,kk_h).*Cu(:,kk_u);
%            sparseWHtUu = sparse(1:noQuadSites,1:noQuadSites,...
%                weightHtUu );
% %            termHtUu = quadDataTensor.BtTr*sparseWHtUu*quadDataTensor.Bu;
%            
%            %Hut'*(...)*Uu
%            weightHutUu = weight5.*Cut(:,kk_h).*Cu(:,kk_u) + ...
%                weight6.*Cuut(:,kk_h).*Cu(:,kk_u) + ...
%                weight7.*Cut(:,kk_h).*Cuu(:,kk_u) + ...
%                weight8.*Cuut(:,kk_h).*Cuu(:,kk_u);
%            sparseWHutUu = sparse(1:noQuadSites,1:noQuadSites,...
%                weightHutUu );
% %            termHutUu = quadDataTensor.ButTr*sparseWHutUu*quadDataTensor.Bu;
%            
%            %Huut'*(...)*Uu
%            weightHuutUu = weight6.*Cut(:,kk_h).*Cu(:,kk_u) + ...
%                weight9.*Cuut(:,kk_h).*Cu(:,kk_u) + ...
%                weight8.*Cut(:,kk_h).*Cuu(:,kk_u);
%            sparseWHuutUu = sparse(1:noQuadSites,1:noQuadSites,...
%                weightHuutUu );
% %            termHuutUu = quadDataTensor.BuutTr*sparseWHuutUu*quadDataTensor.Bu;
%            
%            %Hut'*(...)*Uuu
%            weightHutUuu = weight7.*Cut(:,kk_h).*Cu(:,kk_u) + ...
%                weight8.*Cuut(:,kk_h).*Cu(:,kk_u);
%            sparseWHutUuu = sparse(1:noQuadSites,1:noQuadSites,...
%                weightHutUuu );
% %            termHutUuu = quadDataTensor.ButTr*sparseWHutUuu*quadDataTensor.Buu;
%            
%            %Huut'*(...)*Uuu
%            weightHuutUuu = weight8.*Cut(:,kk_h).*Cu(:,kk_u);
%            sparseWHuutUuu = sparse(1:noQuadSites,1:noQuadSites,...
%                weightHuutUuu );
% %            termHuutUuu = quadDataTensor.BuutTr*sparseWHuutUuu*quadDataTensor.Buu;
%                
%            totalterms = (1/2*quadDataTensor.BuTr*sparseWHuUu + ...
%                quadDataTensor.BuuTr*sparseWHuuUu + ...
%                quadDataTensor.BtTr*sparseWHtUu + ...
%                quadDataTensor.ButTr*sparseWHutUu + ...
%                quadDataTensor.BuutTr*sparseWHuutUu)*quadDataTensor.Bu + ...
%                (1/2*quadDataTensor.BuuTr*sparseWHuuUuu + ...
%                quadDataTensor.ButTr*sparseWHutUuu + ...
%                quadDataTensor.BuutTr*sparseWHuutUuu)*quadDataTensor.Buu;
%            
%            Hess( (kk_h-1)*noControlPoints+1:kk_h*noControlPoints,...
%                (kk_u-1)*noControlPoints+1:kk_u*noControlPoints) = totalterms;
% %            temp = termHuuUu + termHtUu + termHutUu + termHuutUu + ...
% %                termHutUuu + termHuutUuu;
% %            Hess( (kk_h-1)*noControlPoints+1:kk_h*noControlPoints,...
% %                (kk_u-1)*noControlPoints+1:kk_u*noControlPoints) ...
% %                = 1/2*(termHuUu + termHuuUuu) + temp;
%        end
%     end
%     Hess = Hess + Hess';
%     
%     %Symmetric
%     sparseWHuUu = sparse(1:noQuadSites,1:noQuadSites, term1 );
%     termHuUu = quadDataTensor.BuTr*sparseWHuUu*quadDataTensor.Bu;
%         
%     sparseWHtUt = sparse(1:noQuadSites,1:noQuadSites, term3 );
%     termHtUt = quadDataTensor.BtTr*sparseWHtUt*quadDataTensor.Bt;
%     
%     sparseWHutUut = sparse(1:noQuadSites,1:noQuadSites, term4 );
%     termHutUut = quadDataTensor.ButTr*sparseWHutUut*quadDataTensor.But;
%     
%     sparseWHuutUuut = sparse(1:noQuadSites,1:noQuadSites, term6 );
%     termHuutUuut = quadDataTensor.BuutTr*sparseWHuutUuut*quadDataTensor.Buut;
%     
%     %Mixed
%     sparseWHuuUu = sparse(1:noQuadSites,1:noQuadSites, term2 );
%     termHuuUu = quadDataTensor.BuuTr*sparseWHuuUu*quadDataTensor.Bu;
%     
%     sparseWHuutUut = sparse(1:noQuadSites,1:noQuadSites, term5 );
%     termHuutUut = quadDataTensor.BuutTr*sparseWHuutUut*quadDataTensor.But;
%     
%     temp = termHuuUu + termHuutUut;
%     temp2 = termHuUu + termHtUt + termHutUut + termHuutUuut + temp + temp';
%     for kk = 1:splineData.dSpace
%         Hess2( (kk-1)*noControlPoints+1:kk*noControlPoints,...
%                (kk-1)*noControlPoints+1:kk*noControlPoints) ...
%                = temp2;
%     end
%     Hess = Hess + Hess2;
%     
%     
%     %Compute Hessian terms w.r.t. phi, translation, rotation and shift
%     c1Hess = zeros(Nphi+dSpace+2,...
%             Nphi+dSpace+2, ...
%             N*dSpace);    
% %
% %     d1PhiInterpolMat = sparse( [1:noInterpolPoints,(1:noInterpolPoints)+noInterpolPoints],...
% %         repmat(1:noInterpolPoints,1,dSpace),d1PhiInterpol(:)  );
% %     d1PhiDotPhiInterpol = d1PhiInterpolMat*quadData.B_interpolPhi;
% %     d1PhiDotPhiInterpol = reshape( d1PhiDotPhiInterpol, [], dSpace*splineData.Nphi );
% %     %     temp1 = [temp, - d1_0];
%     
%     %c1 dphidphi
%     if optDiff
%         d1uuPhiInterpolMat = sparse( [1:noInterpolPoints,(1:noInterpolPoints)+noInterpolPoints],...
%         repmat(1:noInterpolPoints,1,dSpace),d1uuPhiInterpol(:)  );
%         d1uuPhiDotPhiInterpol = d1uuPhiInterpolMat*quadData.B_interpolPhi;
%         
%         for ii = 1:splineData.Nphi %TODO: Change for-loop to something faster
%             for jj = ii:splineData.Nphi
% %                 Phi_jj_Mat = sparse(1:dSpace*noInterpolPoints,1:dSpace*noInterpolPoints,...
% %                     repmat(quadData.B_interpolPhi(:,jj),dSpace,1));
% %                 d1uuPhiDotPhiPhiInterpol = Phi_jj_Mat*d1uuPhiDotPhiInterpol;
% %                 d1uuPhiDotPhiPhiInterpol = reshape( d1uuPhiDotPhiPhiInterpol, [], dSpace*splineData.Nphi );
% %                 d1uuPhiDotPhiPhi = quadData.B_interpolS \ d1uuPhiDotPhiPhiInterpol;
% %                 d1uuPhiDotPhiPhi = reshape( d1uuPhiDotPhiPhi, splineData.dSpace*splineData.N, [] ); %|[]|=Nphi
% %                 
% %                 c1Hess(ii,jj,:) = d1uuPhiDotPhiPhi(:,jj);
%                 
%                 BijTemp = quadData.B_interpolPhi(:,ii).*quadData.B_interpolPhi(:,jj);
%                 BijTempMat = sparse(1:noInterpolPoints,1:noInterpolPoints,BijTemp);
%                 d1uuPhiDotPhiPhiInterpolTemp = BijTempMat*d1uuPhiInterpol;
%                 d1uuPhiDotPhiPhiTemp = quadData.B_interpolS \ d1uuPhiDotPhiPhiInterpolTemp;
%                 d1uuPhiDotPhiPhiTemp = reshape( d1uuPhiDotPhiPhiTemp, dSpace*N, [] ); %|[]|=Nphi
%                 
%                 c1Hess(ii,jj,:) = d1uuPhiDotPhiPhiTemp;
%                 
%                 c1Hess(jj,ii,:) = c1Hess(ii,jj,:);
%             end
%         end
%         
%         %c1 dphidalpha and dphidbeta
%         d1uuPhiDotPhiInterpol = reshape( d1uuPhiDotPhiInterpol, [], dSpace*Nphi );
%         d1uuPhiDotPhi = quadData.B_interpolS \ d1uuPhiDotPhiInterpol;
%         d1uuPhiDotPhi = reshape( d1uuPhiDotPhi, splineData.dSpace*splineData.N, [] );
%         
%         d1uPhiRot90Interpol = d1uPhiInterpol*rotation90;
%         d1uPhiRot90InterpolMat = sparse( [1:noInterpolPoints,(1:noInterpolPoints)+noInterpolPoints],...
%             repmat(1:noInterpolPoints,1,dSpace),d1uPhiRot90Interpol(:)  );
%         d1uPhiRot90DotPhiInterpol = d1uPhiRot90InterpolMat*quadData.B_interpolPhi;
%         d1uPhiRot90DotPhiInterpol = reshape( d1uPhiRot90DotPhiInterpol, [], dSpace*Nphi );
%         d1uPhiRot90DotPhi = quadData.B_interpolS \ d1uPhiRot90DotPhiInterpol;
%         d1uPhiRot90DotPhi = reshape( d1uPhiRot90DotPhi, dSpace*N, [] ); %|[]|=Nphi
%         
%         if optShift
%             for ii = 1:Nphi
%                 %dphidalpha
%                 c1Hess(ii,end,:) = -d1uuPhiDotPhi(:,ii);
%                 c1Hess(end,ii,:) = c1Hess(ii,end,:);
%             end
%         end
%         if optRot
%             for ii = 1:Nphi
%                 %dphidbeta
%                 c1Hess(ii,end-1,:) = d1uPhiRot90DotPhi(:,ii);
%                 c1Hess(end-1,ii,:) = c1Hess(ii,end-1,:);
%             end
%         end
%     end %dphi terms
%     
%     
%     %c1 dbetadalpha
%     if optRot && optShift
%         c1Hess(end,end-1,:) = reshape(d1uPhiDotAlpha*rotation90,[],1);
%         c1Hess(end-1,end,:) = c1Hess(end,end-1,:);
%     end
%     
%     %c1 dbetadbeta
%     if optRot
% %         d1Phi = quadData.B_interpolS \ d1PhiInterpol;
% %         d1Phi = reshape( d1Phi, splineData.dSpace*splineData.N, [] ); %|[]|=Nphi
%         c1Hess(end-1,end-1,:) = -reshape(d1,dSpace*N,[]);
%     end
%     
%     %c1 dalphadalpha
%     if optShift
%         d1uuPhi = quadData.B_interpolS \ d1uuPhiInterpol;
%         c1Hess(end,end,:) = d1uuPhi(:);
%     end
%     
%     %c1 dvdbeta
%     if optRot && optTra 
%         c1Hess(end-3,end-1,:) = reshape([ones(N,1),zeros(N,1)]*rotation90*rotation,...
%             splineData.dSpace*splineData.N,[]);
%         c1Hess(end-2,end-1,:) = reshape([zeros(N,1),ones(N,1)]*rotation90*rotation,...
%             splineData.dSpace*splineData.N,[]);
%         
%         
%         c1Hess(end-1,end-3,:) = c1Hess(end-3,end-1,:);
%         c1Hess(end-1,end-2,:) = c1Hess(end-2,end-1,:);
%     end
%     
%     %Extract indices for d1 controls in H
%     %Harcoded for dSpace = 2
%     %TODO: For any dSpace.
%     d1Indices  = (noControlPoints-N+1:noControlPoints)';
%     d1Indices = [d1Indices,d1Indices+noControlPoints];
%     d1Indices = d1Indices(:);
%     
%     %Extract dEdc1 parts of control points Hessian
%     dEc1Hess = Hess(d1Indices,d1Indices);
%     
%     %Compute dEdGamma part of Hessian
%     HessGamma = c1dGamma'*dEc1Hess*c1dGamma; %+ ...
%     
% %     Edc1Mat = sparse(1:noInterpolPoints,1:noInterpolPoints,Edc1);
%     c1Hess = permute( c1Hess, [3 1 2]);
%     HessGamma2 = zeros(1,noGamma,noGamma);
%     for ii = 1:noGamma;
%         HessGamma2(:,:,ii) = Edc1*c1Hess(:,:,ii);
%     end
%     HessGamma2 = squeeze(HessGamma2); %Remove singleton dimension
%     
% %     c1Hess = permute( c1Hess,[2 3 1]) ;
% %     HessGamma3 = zeros(noGamma,noGamma);
% %     for ii = 1:noGamma
% %         for jj = ii:noGamma
% %             HessGamma3(ii,jj) = Edc1*squeeze(c1Hess(ii,jj,:));
% %             HessGamma3(jj,ii) = HessGamma3(ii,jj);
% %         end
% %     end
%     
%     %Extract indices for interior controls in H, upper left block
%     intIndices  = zeros( 1,noVariables );
%     for kk = 1:splineData.dSpace
%         intIndices( (kk-1)*noVariables+1:kk*noVariables) = ...
%             (kk-1)*noControlPoints+[splineData.N+1:noControlPoints-splineData.N];
%     end
%         
%     HessInt = Hess( intIndices, intIndices);
%     
%     HessIntd1 = Hess( intIndices, d1Indices);
%     
%     HessIntGamma = HessIntd1*c1dGamma;
%     
%     HessFull = [HessInt,HessIntGamma;HessIntGamma',HessGamma + HessGamma2];
% 
%     H = HessFull;
%     
% %     E = d1;
% %     dE = dEc1Hess;
% %     dE = c1dGamma;
% %     H = permute( c1Hess,[2 3 1]) ;
% %     H = HessGamma + HessGamma2;
% %     H = dEc1Hess;
% %     H = Hess;
% 
% end



end

