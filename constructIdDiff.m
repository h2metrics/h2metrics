%% constructIdDiff
% Constructs identity diffeomoprhism for given splineData. In fact it
% constructs the control points of the periodic function f(x), such that
%   phi(x) = x + f(x)
% is the required diffeomorphism.
%
% Input
%   splineData
%       Determines Nphi
%
% Output
%   phi
%       Control points of the function f.
function [ phi ] = constructIdDiff( splineData )

phi = zeros([splineData.Nphi, 1]);

end

