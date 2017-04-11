function [E, dE] = varifoldDistanceSquared(u, v, splineData)

varData = splineData.varData;
B_varS = splineData.quadData.B_varS;

ptsU = B_varS * u;
ptsV = B_varS * v;

varU = struct( 'x', ptsU, 'G', varData.G );
varV = struct( 'x', ptsV, 'G', varData.G );

objfun.kernel_geom = varData.kernelGeom;
objfun.kernel_grass = varData.kernelGrass;
objfun.kernel_size_geom = varData.kernelSize;

E = varifoldnorm(varU, varV, objfun);

if nargout < 2
    return
end

%% Compute gradient using the chain rule
gradPts = dvarifoldnorm(varU, varV, objfun);
dE = B_varS' * gradPts;

end