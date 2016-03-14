%% pathFlatH1H2Dist
%
% Computes the norm \| d \|_{H1H2} of a spline path.
%
% Input
%   d
%       The spline path
%   splineData
%       General information about the splines used.
%
% Output
%   dist
%       The norm of d
%   comp
%       The six components of the norm separately
%           [ L2L2, L2H1, L2H2, H1L2, H1H1, H1H2 ]
%       We have the identity
%           dist = sqrt( sum( comp.^2 ) )
%
function [ dist, comp ] = pathFlatH1H2Norm( d, splineData, varargin )
    
    %% Get the quadrature data of the spline
    quadData = splineData.quadData;
    
    quadPointsS = quadData.quadPointsS;
    quadWeightsS = quadData.quadWeightsS;
    
    quadPointsT = quadData.quadPointsT;
    quadWeightsT = quadData.quadWeightsT;
 
    quadWeights = reshape(quadWeightsS*quadWeightsT',[],1);
    
    %% Set constants
    a = [1 1 1;1 1 1];
%     if ~isempty(splineData.a)
%         a = [splineData.a; splineData.a];
%     end
    
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
    
    %% Evaluate First spline
    
    B = createTensorCollocationMatrix(quadPointsS, quadPointsT,1,1,splineData);
    Bu = createTensorCollocationMatrix(quadPointsS, quadPointsT,2,1,splineData);
    Buu = createTensorCollocationMatrix(quadPointsS, quadPointsT,3,1,splineData);
    Bt = createTensorCollocationMatrix(quadPointsS, quadPointsT,1,2,splineData);
    But = createTensorCollocationMatrix(quadPointsS, quadPointsT,2,2,splineData);
    Buut = createTensorCollocationMatrix(quadPointsS, quadPointsT,3,2,splineData);
    
    S = B*d;
    Su = Bu*d;
    St = Bt*d;
    Suu = Buu*d;
    Sut = But*d;
    Suut = Buut*d;
    
    %% Evaluate distance
    
    L2L2 = sqrt(quadWeights'*sum( (S).^2,2));
    L2H1 = sqrt(quadWeights'*sum( (Su).^2,2));
    L2H2 = sqrt(quadWeights'*sum( (Suu).^2,2));
    H1L2 = sqrt(quadWeights'*sum( (St).^2,2));
    H1H1 = sqrt(quadWeights'*sum( (Sut).^2,2));
    H1H2 = sqrt(quadWeights'*sum( (Suut).^2,2));
    
    dist = sqrt(a(1,1)*L2L2^2 + a(1,2)*L2H1^2 + a(1,3)*L2H2^2 ...
        + a(2,1)*H1L2^2 + a(2,2)*H1H1^2 + a(2,3)*H1H2^2);
    comp = zeros(6,1);
    if nargout > 1
        comp(1) = a(1,1)*L2L2;
        comp(2) = a(1,2)*L2H1;
        comp(3) = a(1,3)*L2H2;
        comp(4) = a(2,1)*H1L2;
        comp(5) = a(2,2)*H1H1;
        comp(6) = a(2,3)*H1H2;
    end
end

