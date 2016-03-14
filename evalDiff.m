%% evalDiff
%
% Evaluates a diffeomorphism at a given set of points
%
% Input
%   pts
%       Points, were to evaluate the curve
%   phi
%       Control points of the diffeomorphism
%   splineData
%       General information about the splines used.
%
% Output
%   phiPts
%       The diffeomorphism evaluated at given points
%
function phiPts = evalDiff(pts, phi, splineData)

pts = pts(:);
pts = mod(pts, 2*pi);
phiPts = deBoor( splineData.knotsPhi, splineData.nPhi, phi, pts, ...
                    1, 'periodic', true );
phiPts = phiPts + pts;
phiPts = mod(phiPts, 2*pi);

end