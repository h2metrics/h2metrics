%% pathSpline2Spline
%
%
% Computes collocation matrix to evaluate a path or its derivatives
% on a cartesian product of evaluation sites. For interpolation it uses the
% midpoints between knots in the S-direction and Chebyshev points in the
% T-direction.
%
% Input
%   dPath
%       Spline path to interpolate
%   splineData
%       Spline data for dPath
%   splineDataNew
%       Spline data for dPathNew
%
% Output
%   dPathNew
%       Spline representing the same path as dPath with splineDataNew
%
function dPathNew = pathSpline2Spline( dPath, splineData, splineDataNew )

% Choose interpolation sites
innerKnotsS = splineDataNew.innerKnotsS;
if mod(splineDataNew.nS, 2) == 0
    interpolS = innerKnotsS(1:end-1) + 0.5*diff(innerKnotsS);
else
    interpolS = innerKnotsS(1:end-1);
end
interpolT = chbpnt(splineDataNew.knotsT, splineDataNew.nT+1);

% Evaluate dPath at said sites
B = createTensorCollocationMatrix(interpolS, interpolT, 1, 1, splineData);
pts = B * dPath;

% Find control points of new spline
Bnew = createTensorCollocationMatrix(interpolS, interpolT, 1, 1, ...
                                     splineDataNew);
dPathNew = Bnew \ pts;