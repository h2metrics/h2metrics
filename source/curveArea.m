%% curveArea
%
% Calculates the area enclosed by a curve. We use Green's theorem to
% calculate the are. This assumes that the curve is simply connected and
% positively oriented. The formula is
%   area(d) = 1/2 * int_0^{2*pi} x*y' - y*x' d\theta
% This only works for closed curves. Open curves do not enclose an area.
%
% Input
%   d
%       Curve
%   splineData
%       Information about the splines used.
%
% Output
%   A
%       Area
function A = curveArea( d, splineData )

quadData = splineData.quadData;

c = quadData.B_S * d;
c_u = quadData.Bu_S * d;

A = 0.5 * sum( (c(:,1) .* c_u(:,2) - c(:,2) .* c_u(:,1)) .* ...
               quadData.quadWeightsS );

end

