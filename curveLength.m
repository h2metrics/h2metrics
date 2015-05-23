function ell = curveLength(d, splineData, quadData)

cQuad_u = quadData.Bu_S * d;
cSpeed = sum( cQuad_u.^2 , 2).^(1/2);
ell = sum(cSpeed .* quadData.quadWeightsS);

end