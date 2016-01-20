%% curveFlatH2InnerProd
% Computes the norm
%   |u| = \sqrt{ <u,u> },
% where <u,u> is given by curveFlatH2InnerProd.
%
% Input
%   u
%       The spline curve
%   splineData
%       Information about splines used.
%
% Optional inputs
%   a
%       Coefficients for the inner product
%
% Output
%   G
%       The norm
%   components = [L2, H1, H2]
%       We have the identity
%           G = sqrt{L2^2 + H1^2 + H2^2}
%
% Notes
%   The order of precedence for the constants are as follows
%     -) Optional parameter 'a'
%     -) splineData.a
%     -) a = [1 0 1]
%
function [G, components] = curveFlatH2Norm( d, splineData, varargin )

% Set constants to be used in metric
a = [1 0 1];
if ~isempty(splineData.a)
    a = splineData.a;
end

ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'a'
                ii = ii + 1;
                a = varargin{ii};
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1; 
    end
end

quadData = splineData.quadData;

c = quadData.B_S * d;
cu = quadData.Bu_S * d;
cuu = quadData.Buu_S * d;

L2 = sum(sum(c .* c, 2) .* quadData.quadWeightsS);
H1 = sum(sum(cu .* cu, 2) .* quadData.quadWeightsS);
H2 = sum(sum(cuu .* cuu, 2) .* quadData.quadWeightsS);

G = sqrt(a(1) * L2 + a(2) * H1 + a(3) * H2);

if nargout > 1 % Return each component separately
    components = [sqrt(a(1)*L2), sqrt(a(2)*H1), sqrt(a(3)*H2)];
end

end

