function c = evalCurve(pts, d, splineData, deriv)

if nargin < 4
    deriv = 1;
end

c = deBoor( splineData.knotsS, splineData.nS, d, pts, ...
                    deriv, 'periodic', true );
c = c(deriv:deriv:end, :);

end