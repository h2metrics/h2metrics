function [ dist, comp ] = pathFlatH1H2dist( d1,splineData1,quadData1,d2,splineData2,quadData2)
% Computes the H1H2 distance between two spline paths, defined on (possibly)
% two different knot sequences. 
    
    %% Get the quadrature data of the finest spline
    if splineData1.N >= splineData2.N
        quadPointsS = quadData1.quadPointsS;
        quadWeightsS = quadData1.quadWeightsS;
    else
        quadPointsS = quadData2.quadPointsS;
        quadWeightsS = quadData2.quadWeightsS;
    end
    
    if splineData1.Nt >= splineData2.Nt
        quadPointsT = quadData1.quadPointsT;
        quadWeightsT = quadData1.quadWeightsT;
    else
        quadPointsT = quadData2.quadPointsT;
        quadWeightsT = quadData2.quadWeightsT;
    end
    
%     quadPointsS = splineData1.quadPointsS;
%     quadPointsT = splineData1.quadPointsT;
    noQuadPointsS = length(quadPointsS);
    noQuadPointsT = length(quadPointsT);
    
    quadWeights = reshape(quadWeightsS*quadWeightsT',[],1);
    
    %% First spline
    
    B1 = createTensorCollocationMatrix(quadPointsS, quadPointsT,1,1,splineData1);
    B1u = createTensorCollocationMatrix(quadPointsS, quadPointsT,2,1,splineData1);
    B1uu = createTensorCollocationMatrix(quadPointsS, quadPointsT,3,1,splineData1);
    B1t = createTensorCollocationMatrix(quadPointsS, quadPointsT,1,2,splineData1);
    B1ut = createTensorCollocationMatrix(quadPointsS, quadPointsT,2,2,splineData1);
    B1uut = createTensorCollocationMatrix(quadPointsS, quadPointsT,3,2,splineData1);
    
    S1 = B1*d1;
    S1u = B1u*d1;
    S1t = B1t*d1;
    S1uu = B1uu*d1;
    S1ut = B1ut*d1;
    S1uut = B1uut*d1;
    
    %% Second spline
    B2 = createTensorCollocationMatrix(quadPointsS, quadPointsT,1,1,splineData2);
    B2u = createTensorCollocationMatrix(quadPointsS, quadPointsT,2,1,splineData2);
    B2uu = createTensorCollocationMatrix(quadPointsS, quadPointsT,3,1,splineData2);
    B2t = createTensorCollocationMatrix(quadPointsS, quadPointsT,1,2,splineData2);
    B2ut = createTensorCollocationMatrix(quadPointsS, quadPointsT,2,2,splineData2);
    B2uut = createTensorCollocationMatrix(quadPointsS, quadPointsT,3,2,splineData2);
    
    S2 = B2*d2;
    S2u = B2u*d2;
    S2t = B2t*d2;
    S2uu = B2uu*d2;
    S2ut = B2ut*d2;
    S2uut = B2uut*d2;
    
    %% Evaluate distance
    
    L2L2 = sqrt(quadWeights'*sum( (S1 - S2).^2,2));
    L2H1 = sqrt(quadWeights'*sum( (S1u - S2u).^2,2));
    L2H2 = sqrt(quadWeights'*sum( (S1uu - S2uu).^2,2));
    H1L2 = sqrt(quadWeights'*sum( (S1t - S2t).^2,2));
    H1H1 = sqrt(quadWeights'*sum( (S1ut - S2ut).^2,2));
    H1H2 = sqrt(quadWeights'*sum( (S1uut - S2uut).^2,2));
    
    dist = sqrt(L2L2^2 + L2H1^2 + L2H2^2 + H1L2^2 + H1H1^2 + H1H2^2);
    comp = zeros(6,1);
    comp(1) = L2L2;
    comp(2) = L2H1;
    comp(3) = L2H2;
    comp(4) = H1L2;
    comp(5) = H1H1;
    comp(6) = H1H2;
end

