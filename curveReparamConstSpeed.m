function c = curveReparamConstSpeed(d, splineData, quadData)

nS = splineData.nS;
knotsS = splineData.knotsS;
interpolS = splineData.interpolS;
B_interpolS = quadData.B_interpolS;

phiPts = arcLengthInverse(interpolS, d, splineData, quadData);

cPts = deBoor(knotsS, nS, d, phiPts, 1, 'periodic', true);

c = B_interpolS \ cPts;

end