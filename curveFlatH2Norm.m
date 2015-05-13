function [G, components] = curveFlatH2Norm( d, splineData, quadData, ...
                                            varargin)

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

c = quadData.B_S * d;
cu = quadData.Bu_S * d;
cuu = quadData.Buu_S * d;

L2 = sum(sum(c .* c, 2) .* quadData.quadWeightsS);
H1 = sum(sum(cu .* cu, 2) .* quadData.quadWeightsS);
H2 = sum(sum(cuu .* cuu, 2) .* quadData.quadWeightsS);

G = sqrt(a(1) * L2 + a(2) * H1 + a(3) * H2);

if nargout > 1 % Return each component separately
    components = [sqrt(L2), sqrt(H1), sqrt(H2)];
end

end

