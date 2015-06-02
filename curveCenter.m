function [ c, center ] = curveCenter( d, splineData, quadData )

N = splineData.N;
dSpace = splineData.dSpace;

cQuad = quadData.B_S * d;
cQuad_u = quadData.Bu_S * d;
cSpeed = sum( cQuad_u.^2 , 2).^(1/2);
cLength = sum(cSpeed .* quadData.quadWeightsS);

sumLength = 0;
for jj = dSpace:-1:1
    sumLength = sumLength + cLength;
    cCenter(jj, 1) = sum(cQuad(:,jj) .* cSpeed ...
        .* quadData.quadWeightsS) / cLength;
end

c = d - ones([N, 1]) * cCenter';
center = cCenter;

end

