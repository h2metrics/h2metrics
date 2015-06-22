%% energyH2Diff
%
% Computes the energy of dPath after some reparametrizations...
%
% Input
%   dPath
%       Path of curves; matrix of dimension [N*Nt, dSpace]
%   phi
%       Reparametrization of the final curve
%   splineData
%       General information about the splines used.
%   quadData, quadDataTensor
%       Precomputed spline collocation matrices at quadrature points.
%
% Output
%   E
%       Energy of the path.
%
function [H] = energyH2DiffHessian( dPath, phi, v, beta, alpha, ...
    splineData, quadData, quadDataTensor, varargin)

optDiff = true;
optTra = true;
optRot = true;
optShift = true;

%% Optional parameters
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'optdiff'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    optDiff = logical(varargin{ii});
                else
                    error('Invalid value for option ''optDiff''.');
                end
            case 'opttra'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    optTra = logical(varargin{ii});
                else
                    error('Invalid value for option ''optTra''.');
                end
            case 'optrot'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    optRot = logical(varargin{ii});
                else
                    error('Invalid value for option ''optRot''.');
                end
            case 'optshift'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    optShift = logical(varargin{ii});
                else
                    error('Invalid value for option ''optShift''.');
                end
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1; 
    end
end

%% Extract parameters
dSpace = splineData.dSpace;
N = splineData.N;
Nphi = splineData.Nphi;
noGamma = (Nphi + dSpace + 1 + 1); %Assumes dSpace = 2, dim(O(2)) = 1.

a = splineData.a; a0 = a(1); a1 = a(2); a2 = a(3);
noInterpolPoints = size(quadData.B_interpolPhi,1);

rotation = [ cos(beta), sin(beta); ...
                 -sin(beta), cos(beta) ]; %Transpose of rotation matrix
             
noControlPoints = splineData.N*splineData.Nt;
noControls = noControlPoints*splineData.dSpace;
noVariables = splineData.N*(splineData.Nt-2); %BVP
noQuadSites = length( quadDataTensor.quadWeights );

%% Apply diffeomorphism, translation, rotation and shift
d1 = dPath(end-N+1:end,:);

if optDiff && optShift
%     d1 = curveComposeDiff(d1, phi-alpha, splineData, quadData);
    phiPts = quadData.B_interpolPhi * (phi-alpha);
    phiPts = phiPts + splineData.interpolS; %Id  + f
    phiPts = mod(phiPts, 2*pi); %Domain of definition [0,2*pi]
    d1PhiInterpol = deBoorTest(d1,splineData.knotsS, splineData.nS, phiPts,1);
    d1 = quadData.B_interpolS \ d1PhiInterpol;
%     d1 = d1temp;
elseif optShift
    d1 = curveApplyShift(d1, alpha, splineData, quadData);
end
if optRot
    rotation = [ cos(beta), sin(beta); ...
                 -sin(beta), cos(beta) ]; %Transpose of rotation matrix
    d1 = d1 * rotation;
%     d1 = d1temprot;
end
if optTra
    d1 = d1 + ones([N, 1]) * v';
end

dPath2 = dPath;
dPath2(end-N+1:end,:) = d1;

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

%% Remaining stuff (should be identical to energyH2)
%% Energy of the path
% L2 Energy terms
Ct_L2 = CtCt.*Cspeed;

%H1 Energy terms
Ct_H1 = CspeedInv.*CutCut;

%H2 Energy terms
Ct_H2 = CutCut .* CuCuu2 .* CspeedInv7 ...
    - 2 * CutCuut .* CuCuu .* CspeedInv5 ...
    + CuutCuut .* CspeedInv3;

Ct_L2H1H2 = (a(1)*Ct_L2 +a(2)*Ct_H1+ a(3)*Ct_H2);
energyIntegrand = Ct_L2H1H2;

% weightMatrix = spdiags( quadDataTensor.quadWeights );

% Compute final energy
E = quadDataTensor.quadWeights'*energyIntegrand;


%% Compute Gradient
% if nargout > 1
    
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
        
    term1test  = (term1.*Cu(:,kk) + term2.*Cuu(:,kk))'*quadDataTensor.Bu;
    term2test = (term2.*Cu(:,kk))'*quadDataTensor.Buu;
    term3test  = (term3.*Ct(:,kk))'*quadDataTensor.Bt;
    term4test  = (term4.*Cut(:,kk)+term5.*Cuut(:,kk))'*quadDataTensor.But;
    term5test  = (term6.*Cuut(:,kk)+term5.*Cut(:,kk))'*quadDataTensor.Buut;
  
    dE(:,kk) = term1test + term2test + term3test + term4test + term5test;
    end
    
    Edc1 = dE(end-splineData.N+1:end,:);

    dEinner = dE(splineData.N+1:end-splineData.N,:);
    
    %Derivate w.r.t phi, translation, rotation, alpha
    %c1 = exp(i*beta)*d1(phi - alpha) + v
    
    %Translation
    c1dv = [];
    for ii = 1:dSpace
        c1dv = blkdiag(c1dv,ones(N,1));
    end
    
    %Rotation
    rotation90 = [0, 1; -1,0]; %Transpose of rotation matrix
    c1dbeta = (d1 - ones([N, 1]) * v')*rotation90;
    c1dbeta = c1dbeta(:);
    
    %Phi
    
    d1_fixed = dPath(end-N+1:end,:);
    phiPts = quadData.B_interpolPhi * (phi-alpha);
    phiPts = phiPts + splineData.interpolS; %Id  + f
    phiPts = mod(phiPts, 2*pi); %Domain of definition [0,2*pi]
    dNonper = [ d1_fixed; d1_fixed(1:splineData.nS,:) ];    
    
%     B_interpol = spcol( splineData.knotsS, splineData.nS+1, ...
%         brk2knt( phiPts, 3 ),'sparse'); %Up to second derivative
%     B_interpol_per = [ B_interpol(:,1:splineData.nS) + B_interpol(:,end-splineData.nS+1:end), ...
%         B_interpol(:,splineData.nS+1:end-splineData.nS) ];
%     d1_all = B_interpol_per * d1_alt;
%     
    for ii = dSpace:-1:1
        [ knotsS_d,dNonper_d,nS_d] = fastBSplineDer( splineData.knotsS,...
            dNonper(:,ii),splineData.nS);
        [ knotsS_dd,dNonper_dd,nS_dd] = fastBSplineDer( knotsS_d,...
            dNonper_d,nS_d);
        
        d1PhiInterpol(:,ii) = fastBSplineEval(splineData.knotsS, dNonper(:,ii), splineData.nS, phiPts);
        d1uPhiInterpol(:,ii) = fastBSplineEval(knotsS_d, dNonper_d, nS_d, phiPts);
        d1uuPhiInterpol(:,ii) = fastBSplineEval(knotsS_dd, dNonper_dd, nS_dd, phiPts);
    end
    
    %Alternative deBoorTest
    d1All = deBoorTest(d1_fixed,splineData.knotsS, splineData.nS, phiPts,3);
    d1DeBoor = d1All(1:3:end,:);
    d1uDeBoor = d1All(2:3:end,:);
    d1uuDeBoor = d1All(3:3:end,:);
    
    d1PhiInterpol = d1DeBoor;
    d1uPhiInterpol = d1uDeBoor;
    d1uuPhiInterpol = d1uuDeBoor;
    
    d1PhiInterpol = d1PhiInterpol*rotation;
    d1uPhiInterpol = d1uPhiInterpol*rotation;
    d1uuPhiInterpol = d1uuPhiInterpol*rotation;
    
    %d1uPhiDotPhi
    d1uPhiInterpolMat = sparse( [1:noInterpolPoints,(1:noInterpolPoints)+noInterpolPoints],...
        repmat(1:noInterpolPoints,1,dSpace),d1uPhiInterpol(:)  );
    d1uPhiDotPhiInterpol = d1uPhiInterpolMat*quadData.B_interpolPhi;
    d1uPhiDotPhiInterpol = reshape( d1uPhiDotPhiInterpol, [], dSpace*splineData.Nphi );
    %     temp1 = [temp, - d1_0];
    d1uPhiDotPhi = quadData.B_interpolS \ d1uPhiDotPhiInterpol;
    d1uPhiDotPhi = reshape( d1uPhiDotPhi, splineData.dSpace*splineData.N, [] ); %|[]|=Nphi
    c1dphi = d1uPhiDotPhi;
    
    %alpha
    d1uPhiDotAlpha = quadData.B_interpolS \ (-d1uPhiInterpol);
%     d1uPhiDotAlpha = d1uPhiDotAlpha(:);
    c1dalpha = d1uPhiDotAlpha(:);
    
    Edc1 = Edc1(:)';
    
    %TODO: Selection of all toggles for optimiziation
%     if optDiff
        c1dGamma = [ c1dphi,c1dv, c1dbeta, c1dalpha ];
%     else
%         c1dGamma = [ zeros(noInterpolPoints, splineData.Nphi), c1dv, c1dbeta, c1dalpha ];
%     end
    dEdgamma = Edc1*c1dGamma;
    %Collect all dE terms
   
    dE = [dEinner(:); dEdgamma'];
    
% end


%% Compute Hessian
% if nargout > 2
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
    
    
    %Compute Hessian terms w.r.t. phi, translation, rotation and shift
    c1Hess = zeros(splineData.Nphi+splineData.dSpace+2,...
            splineData.Nphi+splineData.dSpace+2, ...
            splineData.N*splineData.dSpace);    
%
%     d1PhiInterpolMat = sparse( [1:noInterpolPoints,(1:noInterpolPoints)+noInterpolPoints],...
%         repmat(1:noInterpolPoints,1,dSpace),d1PhiInterpol(:)  );
%     d1PhiDotPhiInterpol = d1PhiInterpolMat*quadData.B_interpolPhi;
%     d1PhiDotPhiInterpol = reshape( d1PhiDotPhiInterpol, [], dSpace*splineData.Nphi );
%     %     temp1 = [temp, - d1_0];
    

    d1uuPhiInterpolMat = sparse( [1:noInterpolPoints,(1:noInterpolPoints)+noInterpolPoints],...
        repmat(1:noInterpolPoints,1,dSpace),d1uuPhiInterpol(:)  );
    d1uuPhiDotPhiInterpol = d1uuPhiInterpolMat*quadData.B_interpolPhi;
    
    %c1 dphidphi
    for ii = 1:splineData.Nphi %TODO: Change for-loop to something faster
        for jj = ii:splineData.Nphi
%         Phi_jj_Mat = sparse(1:dSpace*noInterpolPoints,1:dSpace*noInterpolPoints,...
%             repmat(quadData.B_interpolPhi(:,jj),dSpace,1));
%         d1uuPhiDotPhiPhiInterpol = Phi_jj_Mat*d1uuPhiDotPhiInterpol;
%         d1uuPhiDotPhiPhiInterpol = reshape( d1uuPhiDotPhiPhiInterpol, [], dSpace*splineData.Nphi );
%         d1uuPhiDotPhiPhi = quadData.B_interpolS \ d1uuPhiDotPhiPhiInterpol;
%         d1uuPhiDotPhiPhi = reshape( d1uuPhiDotPhiPhi, splineData.dSpace*splineData.N, [] ); %|[]|=Nphi
%         
%         c1Hess(ii,jj,:) = d1uuPhiDotPhiPhi(:,jj);
%         
        BijTemp = quadData.B_interpolPhi(:,ii).*quadData.B_interpolPhi(:,jj);
        BijTempMat = sparse(1:noInterpolPoints,1:noInterpolPoints,BijTemp);
        d1uuPhiDotPhiPhiInterpolTemp = BijTempMat*d1uuPhiInterpol;
        d1uuPhiDotPhiPhiTemp = quadData.B_interpolS \ d1uuPhiDotPhiPhiInterpolTemp;
        d1uuPhiDotPhiPhiTemp = reshape( d1uuPhiDotPhiPhiTemp, splineData.dSpace*splineData.N, [] ); %|[]|=Nphi
        
        c1Hess(ii,jj,:) = d1uuPhiDotPhiPhiTemp;
        
        c1Hess(jj,ii,:) = c1Hess(ii,jj,:);
        end
    end
    
    %c1 dphidalpha and dphidbeta
    d1uuPhiDotPhiInterpol = reshape( d1uuPhiDotPhiInterpol, [], dSpace*splineData.Nphi );
    d1uuPhiDotPhi = quadData.B_interpolS \ d1uuPhiDotPhiInterpol;
    d1uuPhiDotPhi = reshape( d1uuPhiDotPhi, splineData.dSpace*splineData.N, [] );
    
    d1uPhiRot90Interpol = d1uPhiInterpol*rotation90;
    d1uPhiRot90InterpolMat = sparse( [1:noInterpolPoints,(1:noInterpolPoints)+noInterpolPoints],...
        repmat(1:noInterpolPoints,1,dSpace),d1uPhiRot90Interpol(:)  );
    d1uPhiRot90DotPhiInterpol = d1uPhiRot90InterpolMat*quadData.B_interpolPhi;
    d1uPhiRot90DotPhiInterpol = reshape( d1uPhiRot90DotPhiInterpol, [], dSpace*splineData.Nphi );
    d1uPhiRot90DotPhi = quadData.B_interpolS \ d1uPhiRot90DotPhiInterpol;
    d1uPhiRot90DotPhi = reshape( d1uPhiRot90DotPhi, splineData.dSpace*splineData.N, [] ); %|[]|=Nphi
    
    for ii = 1:splineData.Nphi
        %dphidalpha
        c1Hess(ii,end,:) = -d1uuPhiDotPhi(:,ii);
        c1Hess(end,ii,:) = c1Hess(ii,end,:);
        
        %dphidbeta
        c1Hess(ii,end-1,:) = d1uPhiRot90DotPhi(:,ii);
        c1Hess(end-1,ii,:) = c1Hess(ii,end-1,:);
    end
    
    %c1 dbetadalpha
    c1Hess(end,end-1,:) = reshape(d1uPhiDotAlpha*rotation90,[],1);
    c1Hess(end-1,end,:) = c1Hess(end,end-1,:);
    
    %c1 dbetadbeta
    d1Phi = quadData.B_interpolS \ d1PhiInterpol;
    d1Phi = reshape( d1Phi, splineData.dSpace*splineData.N, [] ); %|[]|=Nphi
    c1Hess(end-1,end-1,:) = -d1Phi;
    
    %c1 dalphadalpha
    d1uuPhi = quadData.B_interpolS \ d1uuPhiInterpol;
    c1Hess(end,end,:) = d1uuPhi(:);
    
    %Extract indices for d1 controls in H
    %Harcoded for dSpace = 2
    %TODO: For any dSpace.
    d1Indices  = (noControlPoints-N+1:noControlPoints)';
    d1Indices = [d1Indices,d1Indices+noControlPoints];
    d1Indices = d1Indices(:);
    
    dEc1Hess = Hess(d1Indices,d1Indices);
    
    HessGamma = c1dGamma'*dEc1Hess*c1dGamma; %+ ...
    
%     Edc1Mat = sparse(1:noInterpolPoints,1:noInterpolPoints,Edc1);
    
    c1Hess = permute( c1Hess, [3 1 2]);
    HessGamma2 = zeros(1,noGamma,noGamma);
    for ii = 1:noGamma;
        HessGamma2(:,:,ii) = Edc1*c1Hess(:,:,ii);
    end
    HessGamma2 = squeeze(HessGamma2); %Remove singleton dimension
    
%     c1Hess = permute( c1Hess,[2 3 1]) ;
%     HessGamma3 = zeros(noGamma,noGamma);
%     for ii = 1:noGamma
%         for jj = ii:noGamma
%             HessGamma3(ii,jj) = Edc1*squeeze(c1Hess(ii,jj,:));
%             HessGamma3(jj,ii) = HessGamma3(ii,jj);
%         end
%     end
    
    
    %Extract indices for interior controls in H, upper left block
    intIndices  = zeros( 1,noVariables );
    for kk = 1:splineData.dSpace
        intIndices( (kk-1)*noVariables+1:kk*noVariables) = ...
            (kk-1)*noControlPoints+[splineData.N+1:noControlPoints-splineData.N];
    end
        
    HessInt = Hess( intIndices, intIndices);
    
    HessIntd1 = Hess( intIndices, d1Indices);
    
    HessIntGamma = HessIntd1*c1dGamma;
    
    HessFull = [HessInt,HessIntGamma;HessIntGamma',HessGamma + HessGamma2];
%     E = d1;
%     dE = dEc1Hess;
%     dE = c1dGamma;
    H = HessFull;
%     H = permute( c1Hess,[2 3 1]) ;
%     H = HessGamma + HessGamma2;
%     H = dEc1Hess;
%     H = Hess;

% end



end

