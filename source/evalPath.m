%% evalPath
%
% Evaluates a path curves at a given time points.
%
% Input
%   dPath
%       Control points of the path
%   evalT
%       Where to evaluate the path.
%   splineData
%       General information about the splines used.]
%   deriv
%       Which derivative to evaluate; deriv=1 to evaluate curve itself;
%       deriv=2 for first derivative etc.
%
% Output
%   d
%       Control points of the curve dPath(evalT,.); Has dimensions 
%           [N, dSpace, noT]
%
function [ d ] = evalPath(dPath, evalT, splineData, deriv)

if nargin < 4
    deriv = 1;
end

controlPointWeights = spcol( splineData.knotsT, splineData.nT+1, ...
                             brk2knt(evalT, deriv) )';
                         
for jj = splineData.dSpace:-1:1
    d_x = reshape(dPath(:,jj), splineData.N, splineData.Nt) * ...
            controlPointWeights(:, deriv:deriv:end);
    d(:,:,jj) = d_x;
end

d = permute(d, [1 3 2]);

end

