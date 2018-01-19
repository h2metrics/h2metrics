%% geodesicExp
%
% Wrapper function for geodesicForward
%
% Input
%   q0
%       Initial curve
%   v
%       Initial velocity for geodesic
%   splineData
%       General information about the splines used.
%
% Output
%   q
%       Endpoint of discrete geodesic
%
function q = geodesicExp(q0, v, splineData)

q = geodesicForward(q0, q0 + v / splineData.stepsT, splineData.stepsT, ...
                    splineData, 'endpoint');
                
end