function [ G ] = metricMatrix( d, splineData,quadData )
%Compute the metric matrix at the curve d.

a = splineData.a;

%C = quadData.B*Coefs;
Cu = quadData.Bu_S*d;
Cuu = quadData.Buu_S*d;
Cspeed = sum( Cu.^2 , 2).^(1/2);
CspeedInv = 1./Cspeed;
CspeedInv2 = CspeedInv.^2;

CH2First = -( sum(Cu.*Cuu,2) )./Cspeed.^4;
CH2Second = 1./Cspeed.^2;

G = zeros( splineData.N*splineData.dSpace );

for ii = 1:splineData.N
    for jj = ii:splineData.N
        V_i = quadData.B_S(:,ii);
        Vu_i = quadData.Bu_S(:,ii);
        Vuu_i = quadData.Buu_S(:,ii);
        
        V_j = quadData.B_S(:,jj);
        Vu_j = quadData.Bu_S(:,jj);
        Vuu_j = quadData.Buu_S(:,jj);
        
        %L2 Energy terms
        L2 = sum( V_i.*V_j,2);
        
        %H1 Energy terms
        H1 = CspeedInv2.*Vu_i.*Vu_j;
        
        %H2 Energy terms
        V1H2 = CH2First.*Vu_i + CH2Second.*Vuu_i;
        V2H2 = CH2First.*Vu_j + CH2Second.*Vuu_j;
        
        H2 = V1H2.*V2H2;
        
        energyIntegrand = (a(1)*L2+ a(2)*H1 + a(3)*H2).*Cspeed;
        G(ii,jj) = sum(quadData.quadWeightsS.*energyIntegrand);
        G(jj,ii) = G(ii,jj);
    end
end

for kk = 2:splineData.dSpace
    G( (kk-1)*splineData.N+1:kk*splineData.N,...
        (kk-1)*splineData.N+1:kk*splineData.N) = ...
        G(1:splineData.N,1:splineData.N);
end


%Compute final energy
end

% Alternative implementation of metrixMatrix, slower but easier to debug,
% relies on curveRiemH2InnerProd.

% N = splineData.N;
% dSpace = splineData.dSpace;
% 
% G = [];
% for jj = N*dSpace:-1:1
%     for kk = N*dSpace:-1:jj
%         U = zeros([N*dSpace, 1]);
%         V = zeros([N*dSpace, 1]);
%         U(jj) = 1;
%         V(kk) = 1;
%         
%         u = reshape(U, [N, dSpace]);
%         v = reshape(V, [N, dSpace]);
%         
%         G(jj, kk) = curveRiemH2InnerProd(d, u, v, splineData, quadData);
%         G(kk, jj) = G(jj, kk);
%     end
% end

