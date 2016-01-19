%% curveApplyGamma
%
% Transforms the curve by the Gamma structure. The formula is
%   c(th) = R_beta * ( d( phi(th) - alpha ) + v )
% and R_beta = exp(i beta) is the counterclockwise rotation by beta.
%
% Input
%   d
%       Control points of the curve
%   gamma
%       Gamma structures to be applied
%   splineData
%       General information about the splines used.
%   quadData
%       Precomputed spline collocation matrices at quadrature points.
%
% Optional parameters
%   applyDiff = {true (default), false}
%       If false, we ignore gamma.phi
%   applyTra = {true (default), false}
%       If false, we ignore gamma.v
%   applyRot = {true (default), false}
%       If false, we ignore gamma.beta
%   applyShift = {true (default), false}
%       If false, we ignore gamma.alpha
%
% Output
%   c
%       The transformed curve
%   dc
%       Jacobi matrix with respect to non-empty entries of gamma
function [c, dc] = curveApplyGamma(d, gamma, ...
                                   splineData, quadData, varargin)

applyDiff = true;
applyTra = true;
applyRot = true;
applyShift = true;

% Some code for handling optional inputs
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'applydiff'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    applyDiff = logical(varargin{ii});
                else
                    error('Invalid value for option ''applyDiff''.');
                end
            case 'applytra'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    applyTra = logical(varargin{ii});
                else
                    error('Invalid value for option ''applyTra''.');
                end
            case 'applyrot'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    applyRot = logical(varargin{ii});
                else
                    error('Invalid value for option ''applyRot''.');
                end
            case 'applyshift'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    applyShift = logical(varargin{ii});
                else
                    error('Invalid value for option ''applyShift''.');
                end
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1;  
    end
end

applyDiff = applyDiff && ~isempty(gamma.phi);
applyTra = applyTra && ~isempty(gamma.v);
applyRot = applyRot && ~isempty(gamma.beta);
applyShift = applyShift && ~isempty(gamma.alpha);

if applyDiff && applyShift
    phi = gamma.phi - gamma.alpha;
elseif applyDiff
    phi = gamma.phi;
end
if applyShift
    alpha = gamma.alpha;
end
if applyTra
    v = gamma.v;
else
    v = [ 0; 0 ];
end
if applyRot
    rotation = [ cos(gamma.beta), -sin(gamma.beta); ...
                 sin(gamma.beta), cos(gamma.beta) ];
else
    rotation = [ 1, 0; 0, 1 ];
end

if applyDiff
    cPhi = curveComposeDiff(d, phi, splineData, quadData);
elseif applyShift
    cPhi = curveApplyShift(d, alpha, splineData, quadData);
else
    cPhi = d;
end
c = cPhi + ones([splineData.N, 1]) * v';
c = c * rotation';

%% Compute the Jacobi matrix
if nargout == 1
    return
end

N = splineData.N;
Nphi = splineData.Nphi;
dSpace = splineData.dSpace;
interpolS = splineData.interpolS;

if isempty(Nphi)
    Nphi = 1;
end
 
% Translation
dv = zeros(N*dSpace, dSpace);
if applyTra
    for ii = 1:dSpace
        dv((ii-1)*N+1:ii*N,ii) = ones(N,1);
        dv((ii-1)*N+1:ii*N,:) = dv((ii-1)*N+1:ii*N,:) * rotation;
    end
end

% Rotation
dbeta = zeros(N*dSpace, 1);
if applyRot
    rotation90 = [0, -1; 1, 0];
    dbeta = c * rotation90';
    dbeta = reshape(dbeta, [N*dSpace, 1]);
end

% Shift and Diff preparation
if applyDiff
    phiPts = evalDiff(interpolS, phi, splineData);
elseif applyShift
    phiPts = interpolS - alpha;
else
    phiPts = interpolS;
end

phiPts = mod(phiPts, 2*pi);
duPhi = evalCurve(phiPts, d, splineData, 2);
duPhiRot = duPhi * rotation';

% Shift
dalpha = zeros(N*dSpace, 1);
if applyShift
    cuPhiRot = quadData.B_interpolS \ duPhiRot;
    dalpha = reshape(-cuPhiRot, [N*dSpace, 1]);
end

% Diff
dphi = zeros(N*dSpace, Nphi);
if applyDiff
    for ii=1:dSpace
        duPhiRotdPhi = (duPhiRot(:,ii) * ones(1, Nphi)) .* ...
                       quadData.B_interpolPhi;
        dphi((ii-1)*N+1:ii*N,:) = quadData.B_interpolS \ duPhiRotdPhi;
    end
end

%Selection of all toggles for optimiziation
dc = [ logical(applyDiff) * dphi, ...
       logical(applyTra) * dv, ...
       logical(applyRot) * dbeta, ...
       logical(applyShift) * dalpha ];