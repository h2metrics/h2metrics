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
function [Ehor] = energyH2Hor2(dPath,splineData, quadData, quadDataTensor)


%% Extract parameters
dSpace = splineData.dSpace;
N = splineData.N;

a = splineData.a; a0 = a(1); a1 = a(2); a2 = a(3);
noInterpolPoints = size(quadData.B_interpolPhi,1);
    
noControlPoints = splineData.N*splineData.Nt;
noControls = noControlPoints*splineData.dSpace;
noVariables = splineData.N*(splineData.Nt-2); %BVP
noQuadSites = length( quadDataTensor.quadWeights );

%% Apply diffeomorphism, translation, rotation and shift
d1 = dPath(end-N+1:end,:);
dPath2 = dPath;


%% Evaluate path at quadrature sites
Cu = quadDataTensor.Bu*dPath2;
Ct = quadDataTensor.Bt*dPath2;
Cut = quadDataTensor.But*dPath2;
Cuu = quadDataTensor.Buu*dPath2;
Cuuu = quadDataTensor.Buuu*dPath2;
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



%% Energy of the whole path
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
Etotal = quadDataTensor.quadWeights'*energyIntegrand;






B_S = quadData.B_S;
Bu_S = quadData.Bu_S;
Buu_S = quadData.Buu_S;
noQuadS = quadData.noQuadPointsS;
noQuadT = quadData.noQuadPointsT;





%Matrices for vertical projection
V = zeros( noQuadS, splineData.N, splineData.dSpace );
dV = zeros( noQuadS, splineData.N, splineData.dSpace );
ddV = zeros( noQuadS, splineData.N, splineData.dSpace );

Ehor = zeros(noQuadT,1);
A = zeros(splineData.N);
b = zeros(splineData.N,1);
h = zeros( splineData.N, noQuadT );

% arclength derivative terms
H2First = -( sum(Cu.*Cuu,2) )./Cspeed.^4;
H2Second = 1./Cspeed.^2;
H1 = 1./Cspeed;

%Compute pointwise horizontal energy
for tk = 1:noQuadT %Loop over each t_k
    ind = 1+(tk-1)*noQuadS:tk*noQuadS; %Indices for points at t_k
    
    %Compute vertical projection of c_t:
    %P_ver(c_dot) satisfies the weak equations
    %G_c(P_ver(c_dot),v) = G_c(c_dot,v), for all v in Vert
    
    Cdot = Ct(ind,:);
    DsCdot = [H1(ind).*Cut(ind,1), H1(ind).*Cut(ind,2)];
    Ds2Cdot = [H2First(ind).*Cut(ind,1) + H2Second(ind).*Cuut(ind,1),...
        H2First(ind).*Cut(ind,2) + H2Second(ind).*Cuut(ind,2)];
    
    %Create a basis V_i of Vert, and compute spacial derivatives

    
    for i = splineData.N:-1:1
        %V_i = C'*B_i
        V(:,i,1) = Cu(ind,1).*B_S(:,i);
        V(:,i,2) = Cu(ind,2).*B_S(:,i);
        %V_i' = C''*B_i + C'*B_i'
        dV(:,i,1) = Cuu(ind,1).*B_S(:,i) + Cu(ind,1).*Bu_S(:,i);
        dV(:,i,2) = Cuu(ind,2).*B_S(:,i) + Cu(ind,2).*Bu_S(:,i);
        %V_i'' = C'''*B_i + 2*C''*B_i' + C'*B_i''
        ddV(:,i,1) = Cuuu(ind,1).*B_S(:,i) + 2*Cuu(ind,1).*Bu_S(:,i)...
            + Cu(ind,1).*Buu_S(:,i);
        ddV(:,i,2) = Cuuu(ind,2).*B_S(:,i) + 2*Cuu(ind,2).*Bu_S(:,i)...
            + Cu(ind,2).*Buu_S(:,i);
        %1st arclength derivative
        DsV(:,i,1) = H1(ind).*dV(:,i,1);
        DsV(:,i,2) = H1(ind).*dV(:,i,2);
        %2nd arclength derivative
        Ds2V(:,i,1) = H2First(ind).*dV(:,i,1) + H2Second(ind).*ddV(:,i,1);
        Ds2V(:,i,2) = H2First(ind).*dV(:,i,2) + H2Second(ind).*ddV(:,i,2);
    end
    
    %Create linear system
    for i = 1:splineData.N
       for j = i:splineData.N %Symmetric matrix
           %A_ij = G(vi,vj)
           A(i,j) = sum(quadData.quadWeightsS.*Cspeed(ind).*(...
               a(1)*V(:,i,1).*V(:,j,1) + a(1)*V(:,i,2).*V(:,j,2) + a(2)*DsV(:,i,1).*DsV(:,j,1) + a(2)*DsV(:,i,2).*DsV(:,j,2)+...
               a(3)*Ds2V(:,i,1).*Ds2V(:,j,1) + a(3)*Ds2V(:,i,2).*Ds2V(:,j,2)));           
           A(j,i) = A(i,j);
       end
       
       %b(j) = G(cdot,vj);
       b(i) = sum( quadData.quadWeightsS.*Cspeed(ind).*(...
           a(1)*V(:,i,1).*Cdot(:,1) + a(1)*V(:,i,2).*Cdot(:,2) + a(2)*Ds2V(:,i,1).*Ds2Cdot(:,1) + a(2)*Ds2V(:,i,2).*Ds2Cdot(:,2)+ ...
           a(3)*Ds2V(:,i,1).*Ds2Cdot(:,1) + a(3)*Ds2V(:,i,2).*Ds2Cdot(:,2) ));
    end
    
    %Solve for vertical components
    h(:,tk) = A\b;
    
    % Calculate pointwise energy of horizontal component
    %CdotHor = Cdot - [V(:,:,1)*h(:,tk) V(:,:,2)*h(:,tk)];
    %Ds2CdotHor = Ds2Cdot - [Ds2V(:,:,1)*h(:,tk) Ds2V(:,:,2)*h(:,tk)];
    
    %Ehor(tk) = sum( quadData.quadWeightsS.*Cspeed(ind).*(...
    %    sum(CdotHor.*CdotHor,2) + sum(Ds2CdotHor.*Ds2CdotHor,2) ) );
    Evert(tk) = h(:,tk)'*b;
end

%Horizontal energy of the path
Everttotal = sum(quadData.quadWeightsT.*Evert');
Ehor = Etotal-Everttotal;
end

