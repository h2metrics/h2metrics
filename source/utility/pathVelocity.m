%% pathVelocity
%
% Evaluates the derivative of a path at a given time point.
%
% Input
%   dPath
%       Control points of the path
%   evalT
%       Where to evaluate the path.
%   splineData
%       General information about the splines used.
%
% Output
%   d
%       Control points of the curve dPath(evalT,.); Has dimensions 
%           [N, dSpace, noT]
%
function [ d ] = pathVelocity(dPath, evalT, splineData)

d = evalPath(dPath, evalT, splineData, 2);

end

