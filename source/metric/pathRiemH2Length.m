%% pathRiemH2Length
%
% Computes the Riemannian length of a spline path.
%
% Input
%   dPath
%       The spline path
%   splineData
%       General information about the splines used.
%
% Optional inputs
%   a
%       Coefficients for the inner product
%
% Output
%   G
%       The length of the path
%
% Notes
%   The order of precedence for the constants are as follows
%     -) Optional parameter 'a'
%     -) splineData.a
%     -) a = [1 0 1 0 0]
%
function G = pathRiemH2Length( dPath, splineData, varargin )
       
% Set constants to be used in metric
a = [1 0 1 0 0];
if ~isempty(splineData.a)
    a = splineData.a;
end

% Parse optional inputs
p = inputParser;
addParameter(p, 'a', a);
parse(p, varargin{:});
a = p.Results.a;

%% Scale-invariant metrics
scaleInv = 0;
if ~isempty(splineData.scaleInv)
    scaleInv = splineData.scaleInv;
end

if scaleInv
    ell = quadWeights' * Cspeed; % Length of c
    
    % Update coeffecients with the length weights
    a(1) = a(1)/ell^3;
    a(2) = a(2)/ell;
    a(3) = a(3)*ell;
    a(4) = a(4)/ell;
    a(5) = a(5)/ell;
end

%%
quadPtsT = splineData.quadData.quadPointsT;
quadWeightsT = splineData.quadData.quadWeightsT;
noQuadPtsT = splineData.quadData.noQuadPointsT;

C = evalPath(dPath, quadPtsT, splineData);
Ct = pathVelocity(dPath, quadPtsT, splineData);

G = zeros([noQuadPtsT, 1]);
for jj = noQuadPtsT:-1:1
    [G(jj), ~] = ...
        curveRiemH2InnerProd( C(:,:,jj), Ct(:,:,jj), Ct(:,:,jj), ...
                              splineData );
end

G = sum(sqrt(G) .* quadWeightsT);

end