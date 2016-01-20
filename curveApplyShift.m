%% curveApplyShift
%
% Applies a shift by -alpha to the curve. The formula is
%   c(th) = d( th - alpha )
%
% Input
%   d
%       Control points of the curve
%   alpha
%       Shift to be applied
%   splineData
%       General information about the splines used.%
% Output
%   c
%       The transformed curve
function c = curveApplyShift(d, alpha, splineData)

% We don't know, if Nphi is set.
splineData2 = splineData;
if isempty(splineData.Nphi) || isempty(splineData.nPhi)
    splineData2.Nphi = 5;
    splineData2.nPhi = 3;
    splineData2.noInterpolS = 5 * max(splineData2.N, splineData2.Nphi);
    splineData2 = constructKnots(splineData2);
    splineData2 = setupQuadData(splineData2);
end

phi = ones([splineData2.Nphi, 1]) * -alpha;
c = curveComposeDiff( d, phi, splineData2);

end
