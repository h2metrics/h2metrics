function [E, dE1, dE2] = varifoldDistanceSquared(u, v, splineData)

varData = splineData.varData;
B_varS = splineData.quadData.B_varS;

ptsU = B_varS * u;
ptsV = B_varS * v;

varU = struct( 'x', ptsU, 'G', varData.G );
varV = struct( 'x', ptsV, 'G', varData.G );

objfun.kernel_geom = varData.kernelGeom;
objfun.kernel_grass = varData.kernelGrass;
objfun.kernel_size_geom = varData.kernelSizeGeom;
objfun.kernel_size_grass = varData.kernelSizeGrass;


if nargout < 2
    E = varifoldnorm(varU, varV, objfun);
else
    [E, gradPts1, gradPts2] = dvarifoldnorm(varU, varV, objfun);
    dE1 = B_varS' * gradPts1;
    dE2 = B_varS' * gradPts2;
end
    
% elseif nargout == 3   % Compute gradient and Hessian
%     [gradPts , hessPts] = d2varifoldnorm(varU, varV, objfun);
%     dE = B_varS' * gradPts;
% 
%     % Reorder hessian.
%     noPts = splineData.varData.noPts;
%     I1 = 1:2:2*noPts-1;
%     J1 = 1:noPts;
%     K1 = ones(noPts,1);
%     I2 = 2:2:2*noPts;
%     J2 = noPts+1:2*noPts;
%     K2 = ones(noPts,1);
%     
%     B_varSDouble = blkdiag(B_varS, B_varS); 
%     
%     P = sparse([I1,I2], [J1,J2], [K1,K2], 2*noPts, 2*noPts);
%     
%     ddE = B_varSDouble' * P' * hessPts * P * B_varSDouble;

end