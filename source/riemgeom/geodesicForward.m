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

% Do we return only the endpoints?
endpoint = p.Results.endpoint;

% Some useful constants
N = splineData.N;
dSpace = splineData.dSpace;

% Here we save the geodesic
qAll = zeros(N, dSpace, Nsteps+1);
qAll(:,:,1) = q0;
qAll(:,:,2) = q1;

options = optimset('TolFun', 1e-6, ...
                   'Display', 'off');
%options = optimoptions('fsolve');
%options = optimoptions(options,'TolFun', 1e-6);
%options = optimoptions(options,'Display','off');
%options = optimoptions(options,'MaxIter',400);
%options = optimoptions(options,'MaxFunEvals',10000);

for ii = 3:Nsteps+1
    disp(ii);
    q0 = qAll(:,:,ii-2);
    q1 = qAll(:,:,ii-1);
    q2_init = q1 + (q1 - q0);

    F = @(q) LagrangianLeftDer( q, q1, q0, splineData );
    [q2, ~, ~] = fsolve( F, q2_init, options);
    
    qAll(:,:,ii) = q2;
end

if endpoint
    q = qAll(:,:,end);
else
    q = qAll;
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