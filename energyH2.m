function [E, dE, Hess] = energyH2( dPath, splineData, varargin )
%Compute the energy of the path given by coefs and splineData.
%
% Input: dPath, [N*Nt, dSpace], matrix of coefficent of the spline path in
%           vectorized form.
%        splineData, splineData struct
%        quadData, quadData struct (see setupQuadData)
%        quadDataTensor, quadDataTensor struct (see setupQuadData)
%
% Output: E, the computed energy
%         dE, (optional) the gradient of E w.r.t to the inner control points
%         Hess, (optional) the Hessian of E w.r.t to the inner control points

endpointsFixed = 'both';

%% Optional inputs
ii = 1;
while ii <= length(varargin)
    if strcmpi(varargin{ii}, 'endpointsfixed')
        ii = ii + 1;
        endpointsFixed = varargin{ii};
    end
    ii = ii + 1;  
end

%Number of different point types
% noControlPoints = size( dPath,1);

N = splineData.N;
Nt = splineData.Nt;
dSpace = splineData.dSpace;
quadDataTensor = splineData.quadDataTensor;

noControlPoints = N * Nt;
noControls = noControlPoints * dSpace;
noQuadSites = length( quadDataTensor.quadWeights );

a = splineData.a;
a0 = a(1); a1 = a(2); a2 = a(3);

%% Preliminary calculations
Cu = quadDataTensor.Bu * dPath;
Ct = quadDataTensor.Bt * dPath;
Cut = quadDataTensor.But * dPath;
Cuu = quadDataTensor.Buu * dPath;
Cuut = quadDataTensor.Buut * dPath;

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
CspeedInv11 = CspeedInv9 .* CspeedInv2;


%% Energy of the path
% L2 Energy terms
Ct_L2 = CtCt .* Cspeed;

%H1 Energy terms
Ct_H1 = CspeedInv .* CutCut;

%H2 Energy terms
Ct_H2 = CutCut .* CuCuu2 .* CspeedInv7 ...
    - 2 * CutCuut .* CuCuu .* CspeedInv5 ...
    + CuutCuut .* CspeedInv3;

Ct_L2H1H2 = (a(1)*Ct_L2 +a(2)*Ct_H1+ a(3)*Ct_H2);
energyIntegrand = Ct_L2H1H2;

% Compute final energy
E = quadDataTensor.quadWeights' * energyIntegrand;

%% Gradient of the energy
if nargout == 1
    return;
end
    
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

switch endpointsFixed
    case 'left'
        dE = dE(N+1:end,:);
    case 'right'
        dE = dE(1:end-N,:);
    case 'both'
        dE = dE(N+1:end-N,:);
end
% Default case is 'none'

%% Compute Hessian-vector product matrix, W = Hess*h
if nargout == 2
    return;
end

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
Hess = spalloc( noControls,noControls, nnzHess*dSpace^2);
Hess2 = spalloc( noControls,noControls, nnzHess*dSpace);

for kk_h = 1:dSpace
   for kk_u = 1:dSpace
       %Symmetric
       %Hu'*(...)*Uu
       weightHuUu = weight1.*Cu(:,kk_h).*Cu(:,kk_u) + ...
           weight2.*Cuu(:,kk_h).*Cu(:,kk_u) + ...
           weight2.*Cu(:,kk_h).*Cuu(:,kk_u) + ...
           weight3.*Cuu(:,kk_h).*Cuu(:,kk_u);
       sparseWHuUu = sparse(1:noQuadSites, 1:noQuadSites, weightHuUu );

       %Huu'*(...)*Uuu
       weightHuuUuu = weight3.*Cu(:,kk_h).*Cu(:,kk_u);
       sparseWHuuUuu = sparse(1:noQuadSites,1:noQuadSites,...
           weightHuuUuu );

       %Mixed
       %Huu'*(...)*Uu
       weightHuuUu = weight2.*Cu(:,kk_h).*Cu(:,kk_u) + ...
           weight3.*Cu(:,kk_h).*Cuu(:,kk_u);
       sparseWHuuUu = sparse(1:noQuadSites,1:noQuadSites,...
           weightHuuUu );

       %Ht'*(...)*Uu
       weightHtUu = weight4.*Ct(:,kk_h).*Cu(:,kk_u);
       sparseWHtUu = sparse(1:noQuadSites,1:noQuadSites,...
           weightHtUu );

       %Hut'*(...)*Uu
       weightHutUu = weight5.*Cut(:,kk_h).*Cu(:,kk_u) + ...
           weight6.*Cuut(:,kk_h).*Cu(:,kk_u) + ...
           weight7.*Cut(:,kk_h).*Cuu(:,kk_u) + ...
           weight8.*Cuut(:,kk_h).*Cuu(:,kk_u);
       sparseWHutUu = sparse(1:noQuadSites,1:noQuadSites,...
           weightHutUu );

       %Huut'*(...)*Uu
       weightHuutUu = weight6.*Cut(:,kk_h).*Cu(:,kk_u) + ...
           weight9.*Cuut(:,kk_h).*Cu(:,kk_u) + ...
           weight8.*Cut(:,kk_h).*Cuu(:,kk_u);
       sparseWHuutUu = sparse(1:noQuadSites,1:noQuadSites,...
           weightHuutUu );

       %Hut'*(...)*Uuu
       weightHutUuu = weight7.*Cut(:,kk_h).*Cu(:,kk_u) + ...
           weight8.*Cuut(:,kk_h).*Cu(:,kk_u);
       sparseWHutUuu = sparse(1:noQuadSites,1:noQuadSites,...
           weightHutUuu );

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

%Extract variable components of H
switch endpointsFixed
    case 'left'
        noVariables = N * (Nt - 1);
        indices = zeros(1, noVariables);
        for kk = 1:dSpace
            indices( (kk-1)*noVariables+1:kk*noVariables ) = ...
                (kk-1)*noControlPoints + [N+1:noControlPoints];
        end
        Hess = Hess(indices, indices);
        
    case 'right'
        noVariables = N * (Nt - 1);
        indices = zeros(1, noVariables);
        for kk = 1:dSpace
            indices( (kk-1)*noVariables+1:kk*noVariables ) = ...
                (kk-1)*noControlPoints + [1:noControlPoints-N];
        end
        Hess = Hess(indices, indices);
        
    case 'both'
        noVariables = N * (Nt - 2);
        indices  = zeros(1, noVariables);
        for kk = 1:dSpace
            indices( (kk-1)*noVariables+1:kk*noVariables ) = ...
                (kk-1)*noControlPoints + [N+1:noControlPoints-N];
        end
        Hess = Hess(indices, indices);  
end
% Default case is 'none'

end

