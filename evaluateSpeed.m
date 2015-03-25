function [ tPathSpeed ] = evaluateSpeed( d, quadData)
%Compute the pointwise energy at the sites given in quadData

Cu = quadData.Bu*d;
Ct = quadData.Bt*d;
Cut = quadData.But*d;
Cuu = quadData.Buu*d;
Cuut = quadData.Buut*d;
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

energyIntegrand2 = reshape( energyIntegrand, length(quadData.weights{1}), length(quadData.weights{2}) );
tPathSpeed = zeros( 1, length(quadData.weights{2}));

for jj = length(quadData.weights{2}):-1:1
    tPathSpeed(jj) = sum(energyIntegrand2(:,jj).*quadData.weights{1}');
end

%Compute the weights in vector form
quadWeightMult = (quadData.weights{2}'*quadData.weights{1})'; %Multiply u t(i)
weightMatrix = spdiags( quadWeightMult(:), 0, numel(quadWeightMult), numel(quadWeightMult));

%Compute final energy by matrix mulitplication, 
%TODO: use loop instead of sparsematrix mult
E = sum(weightMatrix*energyIntegrand);

end

