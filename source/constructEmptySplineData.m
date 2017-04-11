%% constructEmptySplineData
%
% Function constructs an empty splineData struct. Some default parameters
% are being set; others are set to empty arrays.
%
% Default parameters
%   dSpace = 2
%   phiEps = 1e-12
%   a = [1 0 1]
%
% Output
%   splineData
%       The created struct.
%
function [ splineData ] = constructEmptySplineData

splineData = struct('nS', [],...
    'nT', [], ...
    'nPhi', [], ...
    'N', [], ...
    'Nt', [], ...
    'Nphi', [], ...
    'quadDegree', [], ...
    'dSpace', 2, ...
    'phiEps', 1e-12, ...
    'a', [1 0 1], ...
    'knotsS', [], ...
    'knotsT', [], ...
    'knotsPhi', [], ...
    'innerKnotsS', [], ...
    'innerKnotsT', [], ...
    'innerKnotsPhi', [], ...
    'periodic', 1, ...
    'noInterpolS', [], ...
    'interpolS', [], ...
    'stepsT', [], ...
    'quadData', [], ...
    'quadDataTensor', [] );

end
