function [ tPathSpeed ] = evaluateSpeed(task)
%Compute the pointwise energy at the sites given in quadData

Ct      = task.splineData.Bt*task.results.d;
Ctu     = task.splineData.Btu*task.results.d;
Ctuu    = task.splineData.Btuu*task.results.d;
Cu      = task.splineData.Bu*task.results.d;
Cuu     = task.splineData.Buu*task.results.d;
Cspeed  = sum( Cu.^2 , 2).^(1/2);

%Energy of the path

%L2 Energy terms
Ct_L2 = sum( Ct.^2,2);

%H2 Energy terms
CH2First = -( sum(Cu.*Cuu,2) )./Cspeed.^4;
CH2Second = 1./Cspeed.^2;
CtH2X = (CH2First.*Ctu(:,1) + CH2Second.*Ctuu(:,1));
CtH2Y = (CH2First.*Ctu(:,2) + CH2Second.*Ctuu(:,2));
Ct_H2 = CtH2X.^2+CtH2Y.^2;

Ct_L2andH2 = (Ct_L2 + Ct_H2);
energyIntegrand = Ct_L2andH2.*Cspeed;

energyIntegrand2 = reshape( energyIntegrand, length(task.splineData.WU), length(task.splineData.WT) );
tPathSpeed = zeros( 1, length(task.splineData.WT));

for jj = length(task.splineData.WT):-1:1
    tPathSpeed(jj) = sum(energyIntegrand2(:,jj).*task.splineData.WU');
end

%Compute final energy by matrix mulitplication, 
%TODO: use loop instead of sparsematrix mult
% E = sum(task.splineData.W*energyIntegrand);

end

