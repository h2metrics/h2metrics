%% geodesicForward
%
% Computes geodesic forward shooting
%
% Input
%   q0, q1
%       Initial conditions for discrete exponential map
%   Nsteps
%       Number of iterations to compute
%   splineData
%       General information about the splines used.
%
% Optional parameters
%   'endpoint'
%       Return only endpoint of geodesic path
%
% Output
%   q
%       Discrete geodesic path
%
function q = geodesicForward(q0,q1,Nsteps,splineData,varargin)

%TODO: Avoid using "full" command on B.
a0 = splineData.a(1); %L2
a1 = splineData.a(2); %H1
a2 = splineData.a(3); %H2

%Some code for handling optional inputs
ii = 1;
rule = 'left';
endpoint = 0;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
         switch (lower(varargin{ii}))
             case 'left' %Don't return gradient terms from end curves
                 rule = 'left';
             case 'mid'
                 rule = 'mid';
             case 'right'
                 rule = 'right';
             case 'endpoint'
                 endpoint = 1;
             otherwise
                 error('Invalid option: ''%s''.',varargin{ii});
         end     
    end
    ii = ii + 1;  
end

quadData = splineData.quadData;

%Nsteps = 100;
geodesicPoints = zeros(splineData.N,splineData.dSpace,Nsteps);
exitFlags = zeros(1, Nsteps);
geodesicPoints(:,:,1) = q0;
geodesicPoints(:,:,2) = q1;

options_fsolve = optimset('TolFun',1e-6,'Display','off');
%options_fsolve = optimoptions('fsolve');
%options_fsolve = optimoptions(options_fsolve,'TolFun', 1e-6);
%options_fsolve = optimoptions(options_fsolve,'Display','off');
%options_fsolve = optimoptions(options_fsolve,'MaxIter',400);
%options_fsolve = optimoptions(options_fsolve,'MaxFunEvals',10000);

for ii = 3:Nsteps

%             if ~mod(ii/Nsteps*100,10)
%                 disp(['i=',num2str(ii)]);
%             end
        
        q_init = 2*geodesicPoints(:,:,ii-1) - geodesicPoints(:,:,ii-2);
        
        if strcmp(rule,'left')
            F = @(q) LagragianLeftDer(q,geodesicPoints(:,:,ii-1),...
                geodesicPoints(:,:,ii-2),quadData);
        elseif strcmp(rule,'mid')
            F = @(q) LagragianMidDer(q,geodesicPoints(:,:,ii-1),...
                geodesicPoints(:,:,ii-2),quadData);
        elseif strcmp(rule,'right')
            F = @(q) LagragianLeftDer(q,geodesicPoints(:,:,ii-1),...
                geodesicPoints(:,:,ii-2),quadData);
        end
        
%         Fqinit = F(q_init);        
%         disE = @(q)discreteH2Energy( {geodesicPoints(:,:,ii-2),q,q_init},splineData,quadData );        
%         %clear dE_fd dE_nrm;
%         %dE_fd = zeros(14,1);
%         dE_nrm = zeros(14,1);
%         for ll = 1:14
% %             dE_fd(:,:,ll) = centralDiffJac( disE,geodesicPoints(:,:,ii-1),2*10^(-ll));
%             dE_nrm(ll) = norm( Fqinit - centralDiffJac( disE,geodesicPoints(:,:,ii-1),2*10^(-ll)));
%         end
%         semilogy( dE_nrm )
%         
        [geodesicPoints(:,:,ii),fval,exitFlags(ii)] = fsolve( F, q_init, options_fsolve);
end

q = geodesicPoints;
if endpoint
    q = geodesicPoints(:,:,end);
end

function [ Eder ] = LagragianLeftDer( d2,d1,d0,quadData )
%Compute the derivative of the left-discrete energy of the 3 point path 
%(q0,q1,q2) with respect to q1.
% 

%F = D_2(E_disc)
%F(q2,h) =2*( 2*G_q0(q1-q0,h) - 2*G_q1(q2-q1,h) + D_(q1,h)G_*(q2-q1,q2-q1))
%Solve F(q2,h) = 0, for all h = B_i*e_k.
N = size(d0,1);
Eder = zeros( N, 2);
noQuadPointsS = quadData.noQuadPointsS;

%Compute G(d1-d0,d1-d0)
%W(q0,q1) = G_q0(q1-q0,q1-q0)

d_v1 = d1 - d0;
d_v2 = d2 - d1;

B = (quadData.B_S);
Bu = (quadData.Bu_S);
Buu = (quadData.Buu_S);

%C = quadData.B*Coefs;
C0u = Bu*d0;
C0uu = Buu*d0;
C0speed = sum( C0u.^2 , 2).^(1/2);
C0speedInv = 1./C0speed;
C0speedInv2 = C0speedInv.^2;
C0speedInv4 = C0speedInv2.^2;
%C0speedInv6 = C0speedInv.^6;
C0uC0uu = sum(C0u.*C0uu,2);

C1u = Bu*d1;
C1uu = Buu*d1;
C1speed = sum( C1u.^2 , 2).^(1/2);
C1speedInv = 1./C1speed;
C1speedInv2 = C1speedInv.^2;
C1speedInv3 = C1speedInv.*C1speedInv2;
C1speedInv4 = C1speedInv2.^2;
C1speedInv6 = C1speedInv3.^2;
C1uC1uu = sum(C1u.*C1uu,2);

V1 = B*d_v1;
V1u = Bu*d_v1;
V1uu = Buu*d_v1;

V2 = B*d_v2;
V2u = Bu*d_v2;
V2uu = Buu*d_v2;

%L2 Energy terms
L2V2 = sum( V2.^2,2);

%H1 terms, D_S V = 1/|c'|*V'
H1V2 = C1speedInv2.*sum(V2u.^2,2);

V2_C1Ds_X = C1speedInv.*V2u(:,1);
V2_C1Ds_Y = C1speedInv.*V2u(:,2);

%H2 Energy terms
%D_S^2 V = DS2_1*V' + DS2_2*V'', for each curve
C0Ds2_1 = -( C0uC0uu ).*C0speedInv4;
C0Ds2_2 = C0speedInv2;
C1Ds2_1 = -( C1uC1uu ).*C1speedInv4;
C1Ds2_2 = C1speedInv2;

V1_C0Ds2_X = C0Ds2_1.*V1u(:,1) + C0Ds2_2.*V1uu(:,1);
V1_C0Ds2_Y = C0Ds2_1.*V1u(:,2) + C0Ds2_2.*V1uu(:,2);

V2_C1Ds2_X = C1Ds2_1.*V2u(:,1) + C1Ds2_2.*V2uu(:,1);
V2_C1Ds2_Y = C1Ds2_1.*V2u(:,2) + C1Ds2_2.*V2uu(:,2);

H2V2 = V2_C1Ds2_X.^2 + V2_C1Ds2_Y.^2;

%%% G_q0(v1,h)
L2_V1h_1 = sparse(1:noQuadPointsS,1:noQuadPointsS,V1(:,1))*B;
L2_V1h_2 = sparse(1:noQuadPointsS,1:noQuadPointsS,V1(:,2))*B;

H1_V1h_1 = sparse(1:noQuadPointsS,1:noQuadPointsS,C0speedInv2.*V1u(:,1))*Bu;
H1_V1h_2 = sparse(1:noQuadPointsS,1:noQuadPointsS,C0speedInv2.*V1u(:,2))*Bu;

H2_V1h_1 = sparse(1:noQuadPointsS,1:noQuadPointsS,C0Ds2_1.*V1_C0Ds2_X)*Bu + ...
    sparse(1:noQuadPointsS,1:noQuadPointsS,C0Ds2_2.*V1_C0Ds2_X)*Buu;
H2_V1h_2 = sparse(1:noQuadPointsS,1:noQuadPointsS,C0Ds2_1.*V1_C0Ds2_Y)*Bu + ...
    sparse(1:noQuadPointsS,1:noQuadPointsS,C0Ds2_2.*V1_C0Ds2_Y)*Buu;

%%% G_q1(v2,h)
L2_V2h_1 = sparse(1:noQuadPointsS,1:noQuadPointsS,V2(:,1))*B;
L2_V2h_2 = sparse(1:noQuadPointsS,1:noQuadPointsS,V2(:,2))*B;

H1_V2h_1 = sparse(1:noQuadPointsS,1:noQuadPointsS,V2_C1Ds_X.*C1speedInv)*Bu;
H1_V2h_2 = sparse(1:noQuadPointsS,1:noQuadPointsS,V2_C1Ds_Y.*C1speedInv)*Bu;

H2_V2h_1 = sparse(1:noQuadPointsS,1:noQuadPointsS,C1Ds2_1.*V2_C1Ds2_X)*Bu + ...
    sparse(1:noQuadPointsS,1:noQuadPointsS,C1Ds2_2.*V2_C1Ds2_X)*Buu;
H2_V2h_2 = sparse(1:noQuadPointsS,1:noQuadPointsS,C1Ds2_1.*V2_C1Ds2_Y)*Bu + ...
    sparse(1:noQuadPointsS,1:noQuadPointsS,C1Ds2_2.*V2_C1Ds2_Y)*Buu;

%%% D_(q1,h)G(v2,v2)
C1speed_Var_1 = sparse(1:noQuadPointsS,1:noQuadPointsS,C1speedInv.*C1u(:,1))*Bu;
C1speed_Var_2 = sparse(1:noQuadPointsS,1:noQuadPointsS,C1speedInv.*C1u(:,2))*Bu;

V2_C1Ds_Var_X_1 = sparse(1:noQuadPointsS,1:noQuadPointsS,-V2_C1Ds_X.*C1speedInv3.*C1u(:,1).*V2u(:,1))*Bu;
V2_C1Ds_Var_Y_1 = sparse(1:noQuadPointsS,1:noQuadPointsS,-V2_C1Ds_Y.*C1speedInv3.*C1u(:,1).*V2u(:,2))*Bu;
V2_C1Ds_Var_X_2 = sparse(1:noQuadPointsS,1:noQuadPointsS,-V2_C1Ds_X.*C1speedInv3.*C1u(:,2).*V2u(:,1))*Bu;
V2_C1Ds_Var_Y_2 = sparse(1:noQuadPointsS,1:noQuadPointsS,-V2_C1Ds_Y.*C1speedInv3.*C1u(:,2).*V2u(:,2))*Bu;

V2_C1Ds2_Var_X_1 = sparse(1:noQuadPointsS,1:noQuadPointsS,V2_C1Ds2_X.*(4*C1speedInv6.*C1u(:,1).*C1uC1uu.*V2u(:,1)-C1speedInv4.*(C1uu(:,1).*V2u(:,1)+2*C1u(:,1).*V2uu(:,1))) )*Bu...
    - sparse(1:noQuadPointsS,1:noQuadPointsS,V2_C1Ds2_X.*C1speedInv4.*C1u(:,1).*V2u(:,1))*Buu;
V2_C1Ds2_Var_Y_1 = sparse(1:noQuadPointsS,1:noQuadPointsS,V2_C1Ds2_Y.*(4*C1speedInv6.*C1u(:,1).*C1uC1uu.*V2u(:,2)-C1speedInv4.*(C1uu(:,1).*V2u(:,2)+2*C1u(:,1).*V2uu(:,2))) )*Bu...
    - sparse(1:noQuadPointsS,1:noQuadPointsS,V2_C1Ds2_Y.*C1speedInv4.*C1u(:,1).*V2u(:,2))*Buu;
V2_C1Ds2_Var_X_2 = sparse(1:noQuadPointsS,1:noQuadPointsS,V2_C1Ds2_X.*(4*C1speedInv6.*C1u(:,2).*C1uC1uu.*V2u(:,1)-C1speedInv4.*(C1uu(:,2).*V2u(:,1)+2*C1u(:,2).*V2uu(:,1))) )*Bu...
    - sparse(1:noQuadPointsS,1:noQuadPointsS,V2_C1Ds2_X.*C1speedInv4.*C1u(:,2).*V2u(:,1))*Buu;
V2_C1Ds2_Var_Y_2 = sparse(1:noQuadPointsS,1:noQuadPointsS,V2_C1Ds2_Y.*(4*C1speedInv6.*C1u(:,2).*C1uC1uu.*V2u(:,2)-C1speedInv4.*(C1uu(:,2).*V2u(:,2)+2*C1u(:,2).*V2uu(:,2))) )*Bu...
    - sparse(1:noQuadPointsS,1:noQuadPointsS,V2_C1Ds2_Y.*C1speedInv4.*C1u(:,2).*V2u(:,2))*Buu;

    Eder_1_alt = 2*(...
    (2*C0speed.*quadData.quadWeightsS)'*(a0*L2_V1h_1 + a1*H1_V1h_1+ a2*H2_V1h_1)...
    - (2*C1speed.*quadData.quadWeightsS)'*(a0*L2_V2h_1 + a1*H1_V2h_1 + a2*H2_V2h_1)...
    + ((a0*L2V2 +a1*H1V2 + a2*H2V2).*quadData.quadWeightsS)'*C1speed_Var_1...
    + (2*C1speed.*quadData.quadWeightsS)'*( a1*(V2_C1Ds_Var_X_1 + V2_C1Ds_Var_Y_1)...
    + a2*(V2_C1Ds2_Var_X_1+V2_C1Ds2_Var_Y_1)) );

    Eder_2_alt = 2*(...
    (2*C0speed.*quadData.quadWeightsS)'*(a0*L2_V1h_2 + a1*H1_V1h_2+ a2*H2_V1h_2)...
    - (2*C1speed.*quadData.quadWeightsS)'*(a0*L2_V2h_2 + a1*H1_V2h_2 + a2*H2_V2h_2)...
    + ((a0*L2V2 +a1*H1V2 + a2*H2V2).*quadData.quadWeightsS)'*C1speed_Var_2...
    + (2*C1speed.*quadData.quadWeightsS)'*(a1*(V2_C1Ds_Var_X_2 + V2_C1Ds_Var_Y_2)...
    + a2*(V2_C1Ds2_Var_X_2+V2_C1Ds2_Var_Y_2) ) );

    Eder = [Eder_1_alt',Eder_2_alt'];
    
end

function [ Eder ] = LagragianMidDer( d2,d1,d0,quadData )
%Compute the derivative of the discrete energy of the 3 point path (q0,q1,q2)
% with respect to q1.
% 

%E(q2,h) = 2*G_q0(q1-q0,h) + 2*G_q1(q2-q1,h) + D_(q1,h)G_*(q2-q1,q2-q1)
%Solve F(q2,h) = 0, for all h = B_i*e_k.
N = size(d0,1);
Eder = zeros( N, 2);

%Compute G_(d1-d0)/2 (d1-d0,d1-d0)
%W(q0,q1) = G_q0(q1-q0,q1-q0)

d_v1 = d1 - d0;
d_v2 = d2 - d1;

d1mid = (d1+d0)/2;
d2mid = (d2+d1)/2;

B = full(quadData.B_S);
Bu = full(quadData.Bu_S);
Buu = full(quadData.Buu_S);

%C = quadData.B*Coefs;
C0u = Bu*d1mid;
C0uu = Buu*d1mid;
C0speed = sum( C0u.^2 , 2).^(1/2);
C0uC0uu = sum(C0u.*C0uu,2);

C1u = Bu*d2mid;
C1uu = Buu*d2mid;
C1speed = sum( C1u.^2 , 2).^(1/2);
C1uC1uu = sum(C1u.*C1uu,2);

V1 = B*d_v1;
V1u = Bu*d_v1;
V1uu = Buu*d_v1;

V2 = B*d_v2;
V2u = Bu*d_v2;
V2uu = Buu*d_v2;

%L2 Energy terms
L2V1 = sum( V1.^2,2);
L2V2 = sum( V2.^2,2);

%H2 Energy terms
%D_S^2 V = DS2_1*V' + DS2_2*V'', for each curve
C0Ds2_1 = -( sum(C0u.*C0uu,2) )./C0speed.^4;
C0Ds2_2 = 1./C0speed.^2;
C1Ds2_1 = -( sum(C1u.*C1uu,2) )./C1speed.^4;
C1Ds2_2 = 1./C1speed.^2;

V1_C0Ds2_X = C0Ds2_1.*V1u(:,1) + C0Ds2_2.*V1uu(:,1);
V1_C0Ds2_Y = C0Ds2_1.*V1u(:,2) + C0Ds2_2.*V1uu(:,2);

V2_C1Ds2_X = C1Ds2_1.*V2u(:,1) + C1Ds2_2.*V2uu(:,1);
V2_C1Ds2_Y = C1Ds2_1.*V2u(:,2) + C1Ds2_2.*V2uu(:,2);

H2V1 = V1_C0Ds2_X.^2 + V1_C0Ds2_Y.^2;
H2V2 = V2_C1Ds2_X.^2 + V2_C1Ds2_Y.^2;

%
for kk = N:-1:1
    %%% G_q1mid(v1,h)
    L2_V1h(:,kk,1) = V1(:,1).*B(:,kk);
    L2_V1h(:,kk,2) = V1(:,2).*B(:,kk);
    
    h_C0DS2 = C0Ds2_1.*Bu(:,kk) + C0Ds2_2.*Buu(:,kk);
    H2_V1h(:,kk,1) = V1_C0Ds2_X.*h_C0DS2;
    H2_V1h(:,kk,2) = V1_C0Ds2_Y.*h_C0DS2;
    
    %%% G_q2mid(v2,h)
    L2_V2h(:,kk,1) = V2(:,1).*B(:,kk);
    L2_V2h(:,kk,2) = V2(:,2).*B(:,kk);
    
    h_C1DS2 = C1Ds2_1.*Bu(:,kk) + C1Ds2_2.*Buu(:,kk);
    H2_V2h(:,kk,1) = V2_C1Ds2_X.*h_C1DS2;
    H2_V2h(:,kk,2) = V2_C1Ds2_Y.*h_C1DS2; 
    
    %%% D_(q1mid,h/2)G(v1,v1)
    C0speed_Var(:,kk,1) = 1./C0speed.*C0u(:,1).*Bu(:,kk);
    C0speed_Var(:,kk,2) = 1./C0speed.*C0u(:,2).*Bu(:,kk);
    
    V1_C0Ds2_X_Var(:,kk,1) = 4./C0speed.^6.*C0u(:,1).*Bu(:,kk).*C0uC0uu.*V1u(:,1) - ...
        1./C0speed.^4.*( (C0uu(:,1).*Bu(:,kk) + C0u(:,1).*Buu(:,kk)).*V1u(:,1) ...
        + 2*C0u(:,1).*Bu(:,kk).*V1uu(:,1) );
    V1_C0Ds2_Y_Var(:,kk,1) = 4./C0speed.^6.*C0u(:,1).*Bu(:,kk).*C0uC0uu.*V1u(:,2) - ...
        1./C0speed.^4.*( (C0uu(:,1).*Bu(:,kk) + C0u(:,1).*Buu(:,kk)).*V1u(:,2) ...
        + 2*C0u(:,1).*Bu(:,kk).*V1uu(:,2) );
    
    V1_C0Ds2_X_Var(:,kk,2) = 4./C0speed.^6.*C0u(:,2).*Bu(:,kk).*C0uC0uu.*V1u(:,1) - ...
        1./C0speed.^4.*( (C0uu(:,2).*Bu(:,kk) + C0u(:,2).*Buu(:,kk)).*V1u(:,1) ...
        + 2*C0u(:,2).*Bu(:,kk).*V1uu(:,1) );
    V1_C0Ds2_Y_Var(:,kk,2) = 4./C0speed.^6.*C0u(:,2).*Bu(:,kk).*C0uC0uu.*V1u(:,2) - ...
        1./C0speed.^4.*( (C0uu(:,2).*Bu(:,kk) + C0u(:,2).*Buu(:,kk)).*V1u(:,2) ...
        + 2*C0u(:,2).*Bu(:,kk).*V1uu(:,2) );
    
    
    %%% D_(q2mid,h/2)G(v2,v2)
    C1speed_Var(:,kk,1) = 1./C1speed.*C1u(:,1).*Bu(:,kk);
    C1speed_Var(:,kk,2) = 1./C1speed.*C1u(:,2).*Bu(:,kk);
    
    V2_C1Ds2_X_Var(:,kk,1) = 4./C1speed.^6.*C1u(:,1).*Bu(:,kk).*C1uC1uu.*V2u(:,1) - ...
        1./C1speed.^4.*( (C1uu(:,1).*Bu(:,kk) + C1u(:,1).*Buu(:,kk)).*V2u(:,1) ...
        + 2*C1u(:,1).*Bu(:,kk).*V2uu(:,1) );
    V2_C1Ds2_Y_Var(:,kk,1) = 4./C1speed.^6.*C1u(:,1).*Bu(:,kk).*C1uC1uu.*V2u(:,2) - ...
        1./C1speed.^4.*( (C1uu(:,1).*Bu(:,kk) + C1u(:,1).*Buu(:,kk)).*V2u(:,2) ...
        + 2*C1u(:,1).*Bu(:,kk).*V2uu(:,2) );
    
    V2_C1Ds2_X_Var(:,kk,2) = 4./C1speed.^6.*C1u(:,2).*Bu(:,kk).*C1uC1uu.*V2u(:,1) - ...
        1./C1speed.^4.*( (C1uu(:,2).*Bu(:,kk) + C1u(:,2).*Buu(:,kk)).*V2u(:,1) ...
        + 2*C1u(:,2).*Bu(:,kk).*V2uu(:,1) );
    V2_C1Ds2_Y_Var(:,kk,2) = 4./C1speed.^6.*C1u(:,2).*Bu(:,kk).*C1uC1uu.*V2u(:,2) - ...
        1./C1speed.^4.*( (C1uu(:,2).*Bu(:,kk) + C1u(:,2).*Buu(:,kk)).*V2u(:,2) ...
        + 2*C1u(:,2).*Bu(:,kk).*V2uu(:,2) );
    
    %%% Compute F(q2)
    Eder(kk,1) = 2*sum( ( 2*C0speed.*(L2_V1h(:,kk,1) + H2_V1h(:,kk,1)) ...
    - 2*C1speed.*(L2_V2h(:,kk,1) + H2_V2h(:,kk,1)) ...
    + 1/2*C0speed_Var(:,kk,1).*(L2V1 + H2V1) ...
    + C0speed.*( V1_C0Ds2_X.*V1_C0Ds2_X_Var(:,kk,1) + V1_C0Ds2_Y.*V1_C0Ds2_Y_Var(:,kk,1) ) ...
    + 1/2*C1speed_Var(:,kk,1).*(L2V2 + H2V2) ...
    + C1speed.*( V2_C1Ds2_X.*V2_C1Ds2_X_Var(:,kk,1) + V2_C1Ds2_Y.*V2_C1Ds2_Y_Var(:,kk,1) )...
    ).*quadData.quadWeightsS) ;
    
    Eder(kk,2) = 2*sum( ( 2*C0speed.*(L2_V1h(:,kk,2) + H2_V1h(:,kk,2)) ...
    - 2*C1speed.*(L2_V2h(:,kk,2) + H2_V2h(:,kk,2)) ...
    + 1/2*C0speed_Var(:,kk,2).*(L2V1 + H2V1) ...
    + C0speed.*( V1_C0Ds2_X.*V1_C0Ds2_X_Var(:,kk,2) + V1_C0Ds2_Y.*V1_C0Ds2_Y_Var(:,kk,2) ) ...
    + 1/2*C1speed_Var(:,kk,2).*(L2V2 + H2V2) ...
    + C1speed.*( V2_C1Ds2_X.*V2_C1Ds2_X_Var(:,kk,2) + V2_C1Ds2_Y.*V2_C1Ds2_Y_Var(:,kk,2) )...
    ).*quadData.quadWeightsS) ;
end

end




end

