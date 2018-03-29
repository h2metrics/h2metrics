%% geodesicForward
%
% Computes geodesic forward shooting
%
% Input
%   q0, q1
%       Initial conditions for discrete exponential map
%   Nsteps
%       Number of iterations to compute
%   splineData
%       General information about the splines used.
%
% Optional parameters
%   'endpoint'
%       Return only endpoint of geodesic path
%
% Output
%   q
%       Discrete geodesic path
%
function q = geodesicForward(q0, q1, Nsteps, splineData, varargin)

% Handle optional inputs
p = inputParser;
addParameter(p, 'endpoint', 0);
parse(p, varargin{:});

% Do we return only the endpoints
endpoint = p.Results.endpoint;

N = splineData.N;
dSpace = splineData.dSpace;

% Here we save the geodesic
geodesicPoints = zeros(N, dSpace, Nsteps);
geodesicPoints(:,:,1) = q0;
geodesicPoints(:,:,2) = q1;

options = optimset('TolFun', 1e-6, ...
                   'Display', 'off');
%options = optimoptions('fsolve');
%options = optimoptions(options,'TolFun', 1e-6);
%options = optimoptions(options,'Display','off');
%options = optimoptions(options,'MaxIter',400);
%options = optimoptions(options,'MaxFunEvals',10000);

for ii = 3:Nsteps
       
        q_init = 2 * geodesicPoints(:,:,ii-1) - geodesicPoints(:,:,ii-2);
        
        F = @(q) LagrangianLeftDer(q, ...
                                   geodesicPoints(:,:,ii-1),...
                                   geodesicPoints(:,:,ii-2), ...
                                   splineData );
          
        [ geodesicPoints(:,:,ii), ~, ~ ] = fsolve( F, q_init, options);
end

q = geodesicPoints;
if endpoint
    q = geodesicPoints(:,:,end);
end

end

% E(d0, d1, d2) = G_{d0}(d1-d0, d1-d0) + G_{d1}(d2-d1, d2-d2)
% D_{d1} E(...) (h) = 2 * G_{d0}(d1-d0, h) - 2 * G_{d1}(d2-d1, h)
%                      + D_{d1} G_{}(d2-d1,d2-d1) (h)
function Eder = LagrangianLeftDer(d2, d1, d0, splineData)
   [~, ~, dG] = curveRiemH2InnerProd(d1, d2-d1, d2-d1, splineData);
   
   Eder = 2 * curveRiemH2Flat(d0, d1-d0, splineData) ...
            - 2 * curveRiemH2Flat(d1, d2-d1, splineData) ...
            + dG;
end