function [ C ] = evaluateSplinePath( evalS,evalT, d, splineData, varargin )
%Evaluate the tensor product spline given by d and splineData, at the
%sites evalS, evalT.
%
% Output:
%       C, [ length(evalS), 2*length(evalT)]
%        coloumns (i,i+1) are x and y coordinates at evalS(:),evalT(i)
%
%

B = evaluateBasisFunctions( evalS, evalT, splineData, varargin{:});

C = reshape( B*d, length(evalS),2*length(evalT));

end

