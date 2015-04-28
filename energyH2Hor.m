function [ E ] = energyH2Hor( dPath, splineData,quadData,quadDataTensor )
%Compute the energy of the path given by coefs and splineData.
%
% Input: dPath, [N*Nt, dSpace], matrix of coefficent of the spline path in
%           vectorized form.
%        splineData, splineData struct
%        quadData, quadData struct (see setupQuadData)
%        quadDataTensor, quadDataTensor struct (see setupQuadData)
%
% Output: E, the computed energy
%        dE, (optional) the gradient of E w.r.t to the inner control points   

%TODO: Implement a0,a1,a2 weights.

%Some code for handling optional inputs
% i = 1;
% FixedEnds = 0;
% while i <= length(varargin)
%     if (isa(varargin{i},'char'))
%          switch (lower(varargin{i}))
%              case 'fixedends' %Don't return gradient terms from end curves
%                  FixedEnds = 1;
%              case 'periodic'
%                  disp('TODO: Setup periodic conditions')
%              otherwise
%                  error('Invalid option: ''%s''.',varargin{i});
%          end     
%     end
%     i = i + 1;  
% end

%Number of different point types
noQuadS = quadData.noQuadPointsS;
noQuadT = quadData.noQuadPointsT;

% N = splineData.N;
% Nt = splineData.Nt;

%C = quadData.B*Coefs;
Cu = quadDataTensor.Bu*dPath;
Ct = quadDataTensor.Bt*dPath;
Cut = quadDataTensor.But*dPath;
Cuu = quadDataTensor.Buu*dPath;
Cuuu = quadDataTensor.Buu*dPath;
Cuut = quadDataTensor.Buut*dPath;
Cspeed = sum( Cu.^2 , 2).^(1/2);

%Second arclength derivative terms
H2First = -( sum(Cu.*Cuu,2) )./Cspeed.^4;
H2Second = 1./Cspeed.^2;

B_S = quadData.B_S;
Bu_S = quadData.Bu_S;
Buu_S = quadData.Buu_S;

%Matrices for vertical projection
V = zeros( noQuadS, splineData.N, splineData.dSpace );
dV = zeros( noQuadS, splineData.N, splineData.dSpace );
ddV = zeros( noQuadS, splineData.N, splineData.dSpace );

Ehor = zeros( noQuadT,1);
A = zeros(splineData.N);
b = zeros(splineData.N,1);
h = zeros( splineData.N, noQuadT );

%Compute pointwise horizontal energy
for tk = 1:noQuadT %Loop over each t_k
    ind = 1+(tk-1)*noQuadS:tk*noQuadS; %Indices for points at t_k
    
    %Compute vertical projection of c_t:
    %P_ver(c_dot) satisfies the weak equations
    %G_c(P_ver(c_dot),v) = G_c(c_dot,v), for all v in Vert
    
    Cdot = Ct(ind,:);
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
        %2nd arclength derivative
        Ds2V(:,i,1) = H2First(ind).*dV(:,i,1) + H2Second(ind).*ddV(:,i,1);
        Ds2V(:,i,2) = H2First(ind).*dV(:,i,2) + H2Second(ind).*ddV(:,i,2);
    end
    
    %Create linear system
    for i = 1:splineData.N
       for j = i:splineData.N %Symmetric matrix
           %A_ij = G(vi,vj)
           A(i,j) = sum(quadData.quadWeightsS.*Cspeed(ind).*(...
               V(:,i,1).*V(:,j,1) + V(:,i,2).*V(:,j,2) + ...
               Ds2V(:,i,1).*Ds2V(:,j,1) + Ds2V(:,i,2).*Ds2V(:,j,2)));
           
           A(j,i) = A(i,j);
       end
       
       %b(j) = G(cdot,vj);
       b(i) = sum( quadData.quadWeightsS.*Cspeed(ind).*(...
           V(:,i,1).*Cdot(:,1) + V(:,i,2).*Cdot(:,2) + ...
           Ds2V(:,i,1).*Ds2Cdot(:,1) + Ds2V(:,i,2).*Ds2Cdot(:,2) ));
    end
    
    %Solve for vertical components
    h(:,tk) = A\b;
    
    % Calculate pointwise energy of horizontal component
    CdotHor = Cdot - [V(:,:,1)*h(:,tk) V(:,:,2)*h(:,tk)];
    Ds2CdotHor = Ds2Cdot - [Ds2V(:,:,1)*h(:,tk) Ds2V(:,:,2)*h(:,tk)];
    
    Ehor(tk) = sum( quadData.quadWeightsS.*Cspeed(ind).*(...
        sum(CdotHor.*CdotHor,2) + sum(Ds2CdotHor.*Ds2CdotHor,2) ) );
end

%Horizontal energy of the path
E = sum(quadData.quadWeightsT.*Ehor);

end