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
function [E,dE] = energyH2Diff( dPath, phi, v, beta, alpha, ...
    splineData, quadData, quadDataTensor, varargin)

optDiff = false;
optTra = true;
optRot = true;
optShift = false;

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
    d1 = curveComposeDiff(d1, phi-alpha, splineData, quadData);
elseif optShift
    d1 = curveApplyShift(d1, alpha, splineData, quadData);
end

if optRot
    rotation = [ cos(beta), sin(beta); ...
                 -sin(beta), cos(beta) ];
    d1 = d1 * rotation;
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
% Energy of the path
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
        
    term1test  = (term1.*Cu(:,kk) + term2.*Cuu(:,kk))'*quadDataTensor.Bu;
    term2test = (term2.*Cu(:,kk))'*quadDataTensor.Buu;
    term3test  = (term3.*Ct(:,kk))'*quadDataTensor.Bt;
    term4test  = (term4.*Cut(:,kk)+term5.*Cuut(:,kk))'*quadDataTensor.But;
    term5test  = (term6.*Cuut(:,kk)+term5.*Cut(:,kk))'*quadDataTensor.Buut;
  
    dE(:,kk) = term1test + term2test + term3test + term4test + term5test;
    end
    
    Edc1 = dE(end-splineData.N+1:end,:);

    dE = dE(splineData.N+1:end-splineData.N,:);
    
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
    if optDiff
        phiPts = quadData.B_interpolPhi * (phi-alpha);
        phiPts = phiPts + splineData.interpolS; %Id  + f
        phiPts = mod(phiPts, 2*pi); %Domain of definition [0,2*pi]
    else
        phiPts = splineData.interpolS - alpha;
    end
    dNonper = [ d1_fixed; d1_fixed(1:splineData.nS,:) ];    
    
    for ii = dSpace:-1:1
        [ knotsS_d,dNonper_d,nS_d] = fastBSplineDer( splineData.knotsS,...
            dNonper(:,ii),splineData.nS);
        [ knotsS_dd,dNonper_dd,nS_dd] = fastBSplineDer( knotsS_d,...
            dNonper_d,nS_d);
        
%         d1PhiInterpol(:,ii) = fastBSplineEval(splineData.knotsS, dNonper(:,ii), splineData.nS, phiPts);
        d1uPhiInterpol(:,ii) = fastBSplineEval(knotsS_d, dNonper_d, nS_d, phiPts);
%         d1uuPhiInterpol(:,ii) = fastBSplineEval(knotsS_dd, dNonper_dd, nS_dd, phiPts);
    end
    
    %Alternative deBoorTest
%     d1All = deBoorTest(d1_fixed,splineData.knotsS, splineData.nS, phiPts,3);
%     d1DeBoor = d1All(1:3:end,:);
%     d1uDeBoor = d1All(2:3:end,:);
%     d1uuDeBoor = d1All(3:3:end,:);
    
%     d1PhiInterpol = d1DeBoor;
%     d1uPhiInterpol = d1uDeBoor;
%     d1uuPhiInterpol = d1uuDeBoor;
    
%     d1PhiInterpol = d1PhiInterpol*rotation;
    d1uPhiInterpol = d1uPhiInterpol*rotation;
%     d1uuPhiInterpol = d1uuPhiInterpol*rotation;
    
    %d1uPhiDotPhi
    if optDiff
        d1uPhiInterpolMat = sparse( [1:noInterpolPoints,(1:noInterpolPoints)+noInterpolPoints],...
            repmat(1:noInterpolPoints,1,dSpace),d1uPhiInterpol(:)  );
        d1uPhiDotPhiInterpol = d1uPhiInterpolMat*quadData.B_interpolPhi;
        d1uPhiDotPhiInterpol = reshape( d1uPhiDotPhiInterpol, [], dSpace*splineData.Nphi );
        %     temp1 = [temp, - d1_0];
        d1uPhiDotPhi = quadData.B_interpolS \ d1uPhiDotPhiInterpol;
        d1uPhiDotPhi = reshape( d1uPhiDotPhi, splineData.dSpace*splineData.N, [] ); %|[]|=Nphi
        c1dphi = d1uPhiDotPhi;
    else
        c1dphi = zeros( splineData.N*splineData.dSpace, 1);
    end
    
    %alpha
    d1uPhiDotAlpha = quadData.B_interpolS \ (-d1uPhiInterpol);
%     d1uPhiDotAlpha = d1uPhiDotAlpha(:);
    c1dalpha = d1uPhiDotAlpha(:);
    
    Edc1 = Edc1(:)';
    
    %Selection of all toggles for optimiziation
    c1dGamma = [ logical(optDiff)*c1dphi,...
        logical(optTra)*c1dv, ...
        logical(optRot)*c1dbeta, ...
        logical(optShift)*c1dalpha ];

    dEdgamma = Edc1*c1dGamma;
    %Collect all dE terms
    dE = [dE(:); dEdgamma'];
    
end


end

