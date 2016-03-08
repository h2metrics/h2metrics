%% curveSpline2Spline
%
% Computes spline representation of a given spline curve in with respect to
% a different set of splineData. For interpolation it uses midpoints 
% between knots or the knot points depending on the degree of the spline.
%
% Input
%   d
%       Curve to interpolate
%   splineData
%       Spline data for d
%   splineDataNew
%       Spline data for dNew
%
% Output
%   dNew
%       Spline representing the same curve as d with splineDataNew
%
function dNew = curveSpline2Spline( d, splineData, splineDataNew )

% Choose interpolation sites
interpolS = splineDataNew.interpolS;

% Evaluate d at said sites
pts = deBoor( splineData.knotsS, splineData.nS, d, interpolS, 1, ...
              'periodic', true );

% Find control points of new spline
quadDataNew = splineDataNew.quadData;
dNew = quadDataNew.B_interpolS \ pts;

end