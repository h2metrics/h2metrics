%% evalCurve
%
% Evaluates a curve at a given set of points
%
% Input
%   pts
%       Points, were to evaluate the curve
%   d
%       Control points of the curve
%   splineData
%       General information about the splines used.
%   deriv
%       Which derivative to evaluate; deriv=1 to evaluate curve itself;
%       deriv=2 for first derivative etc.
%
% Output
%   c
%       The curve evaluated at given points
%
function c = evalCurve(pts, d, splineData, deriv)

if nargin < 4
    deriv = 1;
end

per = splineData.curveClosed;

c = deBoor( splineData.knotsS, splineData.nS, d, pts, ...
                    deriv, 'periodic', per );
c = c(deriv:deriv:end, :);

end