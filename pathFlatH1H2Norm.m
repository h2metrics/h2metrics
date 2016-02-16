function [ dist, comp ] = pathFlatH1H2Norm( d,splineData)
% Computes the H1H2 norm of a spline path
%
% d:            Coefficient matrix of spline
% splineData:   SplineData for spline
% quadData:     quadrature data for spline

    
    %% Get the quadrature data of the spline
    quadData = splineData.quadData;
    
    quadPointsS = quadData.quadPointsS;
    quadWeightsS = quadData.quadWeightsS;
    
    quadPointsT = quadData.quadPointsT;
    quadWeightsT = quadData.quadWeightsT;
 
    quadWeights = reshape(quadWeightsS*quadWeightsT',[],1);
    
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
    
    dist = sqrt(L2L2^2 + L2H1^2 + L2H2^2 + H1L2^2 + H1H1^2 + H1H2^2);
    comp = zeros(6,1);
    if nargout > 1
        comp(1) = L2L2;
        comp(2) = L2H1;
        comp(3) = L2H2;
        comp(4) = H1L2;
        comp(5) = H1H1;
        comp(6) = H1H2;
    end
end

