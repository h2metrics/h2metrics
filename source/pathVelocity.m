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

controlPointWeights = spcol( splineData.knotsT, splineData.nT+1, ...
                             brk2knt(evalT, 2) )';
                         
for jj = splineData.dSpace:-1:1
    d_x = reshape(dPath(:,jj), splineData.N, splineData.Nt) * ...
            controlPointWeights(:,2);
    d(:,:,jj) = d_x;
end

d = permute(d, [1 3 2]);

end

