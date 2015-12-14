%% curveFlatH2InnerProd
%
% Computes the inner product
%   <u, v> = \int a(1) * u.v + a(2) * u'.v' + a(3) * u''.v'' d\th
%
% Input
%   u, v
%       The spline curves
%   splineData, quadData
%       General supporting files
%
% Optional inputs
%   a
%       Coefficients for the inner product
%
% Output
%   G
%       The inner product
%   components = [L2, H1, H2]
%       We have the identity
%           G = L2 + H1 + H2
%
% Notes
%   The order of precedence for the constants are as follows
%     -) Optional parameter 'a'
%     -) splineData.a
%     -) a = [1 0 1]
%
function [G, components] = curveFlatH2InnerProd( u, v, splineData, ...
                                                 quadData, varargin)

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

U = quadData.B_S * u;
Uu = quadData.Bu_S * u;
Uuu = quadData.Buu_S * u;

V = quadData.B_S * v;
Vu = quadData.Bu_S * v;
Vuu = quadData.Buu_S * v;

L2 = sum(sum(U .* V, 2) .* quadData.quadWeightsS);
H1 = sum(sum(Uu .* Vu, 2) .* quadData.quadWeightsS);
H2 = sum(sum(Uuu .* Vuu, 2) .* quadData.quadWeightsS);

G = a(1) * L2 + a(2) * H1 + a(3) * H2;

if nargout > 1 % Return each component separately
    components = [a(1)*L2, a(2)*H1, a(3)*H2];
end

end