function [E, dE, ddE] = varifoldDistanceSquared(u, v, splineData)

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

if nargout ==  2      % Compute gradient
    gradPts = dvarifoldnorm(varU, varV, objfun);
    dE = B_varS' * gradPts;
elseif nargout == 3   % Compute gradient and Hessian
    [gradPts , hessPts] = d2varifoldnorm(varU, varV, objfun);
    dE = B_varS' * gradPts;

    % Reorder hessian.
    noPts = splineData.varData.noPts;
    I1 = 1:2:2*noPts-1;
    J1 = 1:noPts;
    K1 = ones(noPts,1);
    I2 = 2:2:2*noPts;
    J2 = noPts+1:2*noPts;
    K2 = ones(noPts,1);
    
    B_varSDouble = blkdiag( splineData.quadData.B_varS, splineData.quadData.B_varS); 
    
    P = sparse( [I1,I2],[J1,J2],[K1,K2],2*noPts,2*noPts);
    
    ddE = B_varSDouble'*P'*hessPts*P*B_varSDouble;
end

end