function [ E,dE ] = energyH2VarifoldHANSO( coeff,pars )
% Function for parsing energyH2Varifold to HANSO

[E,dE] = energyH2Varifold( ...
    [pars.d0; reshape(coeff(1:end-3), [],2)],pars.d1,coeff(end-2),coeff(end-1:end),...
    pars.splineData,'optRot',pars.optRot,'optTra',pars.optTra);


end

