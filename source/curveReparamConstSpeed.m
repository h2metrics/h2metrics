%% curveReparamConstSpeed
%
% Computes spline approximation of the constant speed parametrization of a
% given curve
%
% Input
%   d
%       Control points of the curve
%   splineData
%       General information about the splines used.
%
% Output
%   c
%       The constant speed parametrization
%
function c = curveReparamConstSpeed(d, splineData)

interpolS = splineData.interpolS;
B_interpolS = splineData.quadData.B_interpolS;

phiPts = arcLengthInverse(interpolS, d, splineData);
cPts = evalCurve(phiPts, d, splineData);
c = B_interpolS \ cPts;

end

% Computes the inverse of the constant-speed function
function inverseTheta = arcLengthInverse(evalS, d, splineData)

    quadData = splineData.quadData;

    Cu = quadData.Bu_S*d;
    Cspeed = sum(Cu.^2 , 2).^(1/2);

    % Compute arclength of each knot subinterval
    noSubIntervals = length(splineData.innerKnotsS) - 1;
    knotsArclength = sum(reshape( Cspeed .* quadData.quadWeightsS, ...
                                  [], noSubIntervals ), 1);
    knotsArclength = [0, cumsum( knotsArclength )];
    curveLength = knotsArclength(end);
    knotsArclength = 2*pi / curveLength * knotsArclength;

    % Find right subinterval for each point
    knotIndices = discretize(evalS, knotsArclength);

    inverseTheta = zeros(size(evalS));
    for kk = 1:length(evalS)
        knotIndex = knotIndices(kk);
        knotInterval = splineData.innerKnotsS(knotIndex:knotIndex+1);

        % As long as d is a regular curve, fzero guarantees a solution 
        inverseTheta(kk) = fzero(@(x)(arcLength(x) - evalS(kk)), ...
                                 knotInterval);
    end

    % Computes the constant-speed function defined by 
    %   2*pi / length(c) * int_0^t c_speed du
    % Nested function
    function h = arcLength(t) 
        if abs(t - knotInterval(1)) < eps
            h = knotsArclength(knotIndex);
            return
        end

        % Quadrature data for interval [knot, t]
        [quadPts, quadWeights] = ...
            gaussianQuadratureData( [knotInterval(1), t], ...
                                    'degree', splineData.quadDegree(1) );

        % Evaluate curve speed at quadrature sites
        Du = evalCurve(quadPts, d, splineData, 2);
        quadPtsSpeed = sum( Du.^2, 2).^(1/2);

        h = knotsArclength(knotIndex) + ...
                2*pi / curveLength * sum(quadPtsSpeed .* quadWeights');
    end
end
