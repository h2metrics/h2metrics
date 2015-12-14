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
%   applyShift = {true (default), false}
%       If false, we ignore gamma.alpha
%
% Output
%   c
%       The transformed curve
function c = curveApplyGamma(d, gamma, splineData, quadData, varargin)

applyDiff = true;
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

c = d;
if applyDiff && ~isempty(gamma.phi)
    if applyShift && ~isempty(gamma.alpha)
        c = curveComposeDiff( c, gamma.phi - gamma.alpha, ...
                              splineData, quadData );
    else
        c = curveComposeDiff( c, gamma.phi, ...
                              splineData, quadData );
    end
elseif applyShift && ~isempty(gamma.alpha)
    c = curveApplyShift(c, gamma.alpha, splineData, quadData);
end
if ~isempty(gamma.v)
    c = c + ones([splineData.N, 1]) * gamma.v';
end
if ~isempty(gamma.beta)
    rotation = [ cos(gamma.beta), sin(gamma.beta); ...
                 -sin(gamma.beta), cos(gamma.beta) ];
    c = c * rotation;
end
