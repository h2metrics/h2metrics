function [ E, dE ] = energyH2( dPath,splineData,quadDataTensor )
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

%Number of different point types
noControlPoints = size( dPath,1);
N = splineData.N;
%Nt = splineData.Nt;

%C = quadData.B*Coefs;
Cu = quadDataTensor.Bu*dPath;
Ct = quadDataTensor.Bt*dPath;
Cut = quadDataTensor.But*dPath;
Cuu = quadDataTensor.Buu*dPath;
Cuut = quadDataTensor.Buut*dPath;
Cspeed = sum( Cu.^2 , 2).^(1/2);

%Energy of the path

%L2 Energy terms
Ct_L2 = sum( Ct.^2,2);

%H2 Energy terms
CH2First = -( sum(Cu.*Cuu,2) )./Cspeed.^4;
CH2Second = 1./Cspeed.^2;
CtH2X = (CH2First.*Cut(:,1) + CH2Second.*Cuut(:,1));
CtH2Y = (CH2First.*Cut(:,2) + CH2Second.*Cuut(:,2));
Ct_H2 = CtH2X.^2+CtH2Y.^2;

Ct_L2andH2 = (Ct_L2 + Ct_H2);
energyIntegrand = Ct_L2andH2.*Cspeed;

%weightMatrix = spdiags( quadDataTensor.quadWeights );

%Compute final energy
E = sum(quadDataTensor.quadWeights.*energyIntegrand);

if nargout > 1 %Compute gradient  
    FixedEnds = 1;

    % If fixed ends, then ignore first and last N tensor control points
    %controlPointsActive = noControlPoints:-1:1;
    if FixedEnds
        controlPointsActive = noControlPoints-N:-1:N+1;
    end
    
    noActivePoints = length(controlPointsActive);
    noEvalSites = size(quadDataTensor.B,1);
    %The change in the index of a control point in the gradient is
    %hardcoded using N*FixedEnds. This should be changed if other
    %constraints are considered.
    
    %Speed-terms derivative
    gradientSpeedIntegrand = zeros( noEvalSites, 2*noActivePoints);
    for jj = controlPointsActive
        gradientSpeedX = Ct_L2andH2.*(quadDataTensor.Bu(:,jj).*Cu(:,1)./Cspeed);
        gradientSpeedY = Ct_L2andH2.*(quadDataTensor.Bu(:,jj).*Cu(:,2)./Cspeed);
        
        gradientSpeedIntegrand(:,[jj-N*FixedEnds,jj-N*FixedEnds+noActivePoints]) = ...
            [gradientSpeedX, gradientSpeedY];
    end
    
    %L2-terms derivative
    gradientL2Integrand = zeros( noEvalSites, 2*noActivePoints);
    for jj = controlPointsActive
        gradientL2X = 2*quadDataTensor.Bt(:,jj).*Ct(:,1);
        gradientL2Y = 2*quadDataTensor.Bt(:,jj).*Ct(:,2);
        
        gradientL2Integrand(:,[jj-N*FixedEnds,jj-N*FixedEnds+noActivePoints]) = ...
            [gradientL2X.*Cspeed, gradientL2Y.*Cspeed];
    end
    
    
    for jj = controlPointsActive
        %k=1, X-coordinate of control points
        Term1k1 = -2*quadDataTensor.Bu(:,jj).*Cu(:,1)./Cspeed.^4;
        Term2k1 = -(quadDataTensor.Bu(:,jj).*Cuu(:,1)+quadDataTensor.Buu(:,jj).*Cu(:,1))./Cspeed.^4 ...
            + 4*quadDataTensor.Bu(:,jj).*Cu(:,1).*sum(Cu.*Cuu,2)./Cspeed.^6;
        Term3 = (Cspeed.^2.*quadDataTensor.Buut(:,jj) - sum(Cu.*Cuu,2).*quadDataTensor.But(:,jj) )./Cspeed.^4;
        
        Der_CtH2_k1_dijX = Term1k1.*Cuut(:,1) + Term2k1.*Cut(:,1) + Term3;
        Der_CtH2_k1_dijY = Term1k1.*Cuut(:,2) + Term2k1.*Cut(:,2);
        
        gradientH2X = 2*(Der_CtH2_k1_dijX.*CtH2X +  Der_CtH2_k1_dijY.*CtH2Y);
        
        %k=2, Y-coordinate of control points
        Term1k2 = -2*quadDataTensor.Bu(:,jj).*Cu(:,2)./Cspeed.^4;
        Term2k2 = -(quadDataTensor.Bu(:,jj).*Cuu(:,2)+quadDataTensor.Buu(:,jj).*Cu(:,2))./Cspeed.^4 ...
            + 4*quadDataTensor.Bu(:,jj).*Cu(:,2).*sum(Cu.*Cuu,2)./Cspeed.^6;
        
        Der_CtH2_k2_dijX = Term1k2.*Cuut(:,1) + Term2k2.*Cut(:,1);
        Der_CtH2_k2_dijY = Term1k2.*Cuut(:,2) + Term2k2.*Cut(:,2) + Term3;
        gradientH2Y = 2*(Der_CtH2_k2_dijX.*CtH2X +  Der_CtH2_k2_dijY.*CtH2Y);
        
        %Integrand
        gradientH2Integrand(:,[jj-N*FixedEnds,jj-N*FixedEnds+noActivePoints]) = ...
            [gradientH2X.*Cspeed, gradientH2Y.*Cspeed];
        
    end
    
%     keyboard
%     sum( sum( (gradientH2Integrand - gradientH2Integrand_alt).^2 ) )
%     
    for ii = noActivePoints*splineData.dSpace:-1:1
        dE(ii) = sum( quadDataTensor.quadWeights.*(gradientL2Integrand(:,ii) + gradientH2Integrand(:,ii) + gradientSpeedIntegrand(:,ii)) );
%     dE = sum( weightMatrix*(gradientL2H2Integrand + gradientSpeedIntegrand) );
    end
    dE = reshape( dE, noActivePoints, splineData.dSpace);
%     if FixedEnds
%         dE = dE( quadData.noControlPoints(1)+1:end-quadData.noControlPoints(1),:);
%     end
end

end

