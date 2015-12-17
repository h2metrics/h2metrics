%% curveSpline2Spline
%
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
innerKnotsS = splineDataNew.innerKnotsS;
if mod(splineDataNew.nS, 2) == 0 % See deBoor for reasons for this.
    interpolS = innerKnotsS(1:end-1) + 0.5*diff(innerKnotsS);
else
    interpolS = innerKnotsS(1:end-1);
end

% Evaluate d at said sites
pts = deBoor( splineData.knotsS, splineData.nS, d, interpolS, 1, ...
              'periodic', true );

% Find control points of new spline
nS = splineDataNew.nS;
Bnew = spcol( splineDataNew.knotsS, nS+1, ...
              brk2knt( interpolS, 1 ), 'sparse');
Bnew = [ Bnew(:,1:nS) + Bnew(:,end-nS+1:end), ... % For periodic spline
        Bnew(:,nS+1:end-nS) ];
dNew = Bnew \ pts;