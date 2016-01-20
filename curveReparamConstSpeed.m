function c = curveReparamConstSpeed(d, splineData)

quadData = splineData.quadData;
nS = splineData.nS;
knotsS = splineData.knotsS;
interpolS = splineData.interpolS;
B_interpolS = quadData.B_interpolS;

phiPts = arcLengthInverse(interpolS, d, splineData, quadData);

cPts = deBoor(knotsS, nS, d, phiPts, 1, 'periodic', true);

c = B_interpolS \ cPts;

end

function [ inverseTheta ] = arcLengthInverse(evalS, d, splineData, quadData )
% Compute the inverse of the "arc-length function" (the constant speed
% function)

%C = quadData.B*Coefs;
Cu = quadData.Bu_S*d;
%Cuu = quadData.Buu_S*d;
Cspeed = sum( Cu.^2 , 2).^(1/2);

%Compute arclength of each knot subinterval
knotsArclength = sum(reshape(Cspeed.*quadData.quadWeightsS, [], splineData.N),1);
knotsArclength = [0, cumsum( knotsArclength )];
curveLength = knotsArclength(end);
knotsArclength = 2*pi/curveLength*knotsArclength;

inverseTheta = zeros(size(evalS));

%Compute the inverse theta-values
for kk = 1 : length(evalS)
    
    knotIndex = find(knotsArclength == evalS(kk),1);
    if knotIndex %h(kk) is the arclength on a knot.
        inverseTheta(kk) = splineData.innerKnotsS(knotIndex);
        continue %next h-value
    end
    
    %h(i) is not the value on a knot, so as long as the curve is regular, 
    %it is guaranteed to lie in the image of a knot subinterval. fzero is then
    %guaranteed to succeed. 
    knotIndex = find(knotsArclength > evalS(kk),1);
    knotInterval = [splineData.innerKnotsS(knotIndex-1:knotIndex)];
    %As long as d is a regular curve, fzero should guarantee a solution 
    %options = optimset('Display','iter');
    inverseTheta(kk) = fzero(@(x) (arcLength(x) - evalS(kk)),knotInterval);
end


%Nested function
function [h, dh] = arcLength(t) 
%Computes the arclength function (actually the "constant speed on [0,2*pi]" function)
%defined by 2*pi/length(c)*int_t0^t c_speed du
% Input: t, [1 x numTval] , row/col-vector of t-values
%
% Output: h, [1 x numTval], arclength values at t
h = zeros(size(t));

for j = 1 : length(t) %Loop over all t-values
    %Is t(i) a knot? If so, no calculation needed.
    tknot = find( splineData.innerKnotsS == t(j) ); %tknot is an empty matrix if t(i) is not a knot
    if tknot
        h(j) = knotsArclength(tknot);
        continue %Go to next t value
    end
    
    %What knot subinterval is t(i) in? Locate first knot which is larger
    %than t(i), use the one before. 
    knotNo = find(splineData.innerKnotsS > t(j),1) - 1;
    
    %Quadrature data for interval [knot(j) t(i)]
    [quadPts, quadWeights] = gaussianQuadratureData(...
        [splineData.innerKnotsS(knotNo),t(j)],...
        'degree',splineData.quadDegree(1));
    
    %TODO: Use C spline library?
    %Evaluate speed at quadrature sites
    quadPtsEval = spcol( splineData.knotsS, splineData.nS+1, brk2knt( quadPts, 2),'sparse');
    quadPtsEval_per = [quadPtsEval(:,1:splineData.nS) + ...
        quadPtsEval(:,end-splineData.nS+1:end), quadPtsEval(:,splineData.nS+1:end-splineData.nS)];
    quadPtsDer = quadPtsEval_per(2:2:end,:)*d;
    quadPtsSpeed = sum( quadPtsDer.^2 , 2).^(1/2);
    
    h(j) = knotsArclength(knotNo) + 2*pi*sum( quadPtsSpeed.*quadWeights')/curveLength;
end

if nargout > 1
    tEval = spcol( splineData.knotsS, splineData.nS+1, brk2knt( t, 2),'sparse');
    tEval_per = [tEval(:,1:splineData.nS) + ...
        tEval(:,end-splineData.nS+1:end), tEval(:,splineData.nS+1:end-splineData.nS)];
    tDer = tEval_per(2:2:end,:)*d;
    dh = sum( tDer.^2 , 2).^(1/2)/curveLength*2*pi;
end

end % end arclength


end
