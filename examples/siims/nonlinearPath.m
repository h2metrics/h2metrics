%% pathSpline2Spline
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
function dPath = nonlinearPath(d0, d1, splineData)

if ~isequal( size(d0), size(d1) )
    error('Dimension mismatch.');
end

innerKnotsS = splineData.innerKnotsS;
if mod(splineData.nS, 2) == 0
    interpolS = innerKnotsS(1:end-1) + 0.5*diff(innerKnotsS);
else
    interpolS = innerKnotsS(1:end-1);
end
interpolT = chbpnt(splineData.knotsT, splineData.nT+1);
pts=[];
% Evaluate non-linear path at said sites
for jj=1:length(interpolT)
    pts_new =   evalCurve(interpolS,d0, splineData) * (1 -sin(interpolT(jj)*pi/2)) + evalCurve(interpolS,d1, splineData) * (sin(interpolT(jj)*pi/2));
    pts=[pts;pts_new];
end


% Find control points of new spline
B = createTensorCollocationMatrix(interpolS, interpolT, 1, 1, splineData);
dPath = B \ pts;