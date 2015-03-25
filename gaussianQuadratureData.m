function [ Pts, Wts ] = gaussianQuadratureData( Domain, varargin )
%gaussianQuadratureData - compute function arguments and weights for Gaussian quadrature integration.
%
%Usage:
%  [ Pts, Wts ] = gaussianQuadratureData( Domain )
%  [ Pts, Wts ] = gaussianQuadratureData( Domain, 'Degree' , Degree )
%  [ Pts, Wts ] = gaussianQuadratureData( Domain, 'NPoints' , NPoints )
%
%Input:
%  Domain: real, [nIntervals x 1] or [nIntervals x 2].
%     Intervals for the domain of integration. In increasing order 
%
%Options:
%  Degree: integer, [1 x 1].
%    Fixes the quadrature degree (default is 3).
%  NPoints: integer, [1 x 1].
%    Fixes the number of quadrature points (default is 2). 
%
%Output:
%  Pts: real, [(nIntervals x NPoints) x 1] or [1 x (nIntervals x NPoints)].
%    Function arguments.
%  Wts: real, [(nIntervals x NPoints) x 1] or [1 x (nIntervals x NPoints)].
%    Function weights.
%
%Description:
%  Computes function arguments and weights intended for Gaussian
%  quadrature integration. Using these, the Integral of a univariate
%  function Fun(x) of one variable over the domain D = [a,b] may be
%  approximated as Integral ~= Wts' * Fun(Pts), assuming appropriate
%  dimensions of Wts and Fun(Pts). The method is exact for polynomial
%  functions up to Degree.
%
%Reference:
%  * Function arguments and weights are compiled from M. Abromovits &
%    I.A. Stegun (ed.): "Handbook of Mathematical Functions with
%    Formulas, Graphs, and Mathematical Tables", 9th ed.
%
%Note:
%  * Please use the method with caution, as it has not been thoroughly
%    tested yet, and please report bugs to penn@dtu.dk.
%  * The degree and the number of quadrature points fulfill the relation
%    Degree = 2 x NPoints - 1.

% Set/get optional arguments.
NPoints = 2;
Degree  = 2*NPoints - 1;
for i = 1:2:length(varargin)
    if (~ischar(varargin{i}))
        error('Please specify option name as a string.')
    elseif (i+1 > length(varargin))
        error( 'Please specify a value for this option: %s', varargin{i} )
    end
    if (strcmpi(varargin{i}, 'degree'))
        Degree  = varargin{i+1};
        NPoints = ceil((Degree+1)/2); % Adjust the number of points
    elseif (strcmpi(varargin{i}, 'npoints'))
        NPoints = varargin{i+1};
    else
        error('No such option: %s', varargin{i} )
    end
end

% Make a quick and dirty validation of the input.
if (~isnumeric(Domain) || numel(Domain)<2 || any(diff(Domain)<0) )
    error('Please specify interval correctly.');
elseif (~isnumeric(NPoints) || numel(NPoints)~=1)
    error('Please specify the number of points correctly.');
elseif (~isnumeric(Degree) || numel(Degree)~=1)
    error('Please specify degree correctly.');
end

% Get the arguments and weights on the interval [-1,1].
DIM        = size(Domain);
DIM(DIM>1) = NPoints; % Ensure similar dimension as input: [a,b] vs. [a,b]'
nodes      = nan(DIM);
weights    = nan(DIM);

switch (NPoints)
    case (1)
        nodes       = 0.;
        weights     = 2.;
    case (2)
        ii          = [1,2];
        nodes(ii)   = 0.577350269189626*[-1. 1.];
        weights(ii) = 1.               *[1.  1.];
    case (3)
        ii          = [2,3,1];
        nodes(ii)   = [0.                         0.774596669241483*[1. -1.]];
        weights(ii) = [0.888888888888889          0.555555555555556*[1.  1.]];
    case (4)
        ii          = [3,2,4,1];
        nodes(ii)   = [0.339981043584856*[1. -1.] 0.861136311594053*[1. -1.]];
        weights(ii) = [0.652145154862546*[1.  1.] 0.347854845137454*[1.  1.]];
    case (5)
        ii          = [3,4,2,5,1];
        nodes(ii)   = [0                          0.538469310105683*[1. -1.] 0.906179845938664*[1. -1.]];
        weights(ii) = [0.568888888888889          0.478628670499366*[1.  1.] 0.236926885056189*[1.  1.]];
    case (6)
        ii          = [4,3,5,2,6,1];
        nodes(ii)   = [0.238619186083197*[1. -1.] 0.661209386466265*[1. -1.] 0.932469514203152*[1. -1.]];
        weights(ii) = [0.467913934572691*[1.  1.] 0.360761573048139*[1.  1.] 0.171324492379170*[1.  1.]];
    case (7)
        ii          = [4,5,3,6,2,7,1];
        nodes(ii)   = [0.                         0.405845151377397*[1. -1.] 0.741531185599394*[1. -1.] 0.949107912342759*[1. -1.]];
        weights(ii) = [0.417959183673469          0.381830050505119*[1.  1.] 0.279705391489277*[1.  1.] 0.129484966168870*[1.  1.]];
    case (8)
        ii          = [5,4,6,3,7,2,8,1];
        nodes(ii)   = [0.183434642495650*[1. -1.] 0.525532409916329*[1. -1.] 0.796666477413627*[1. -1.] 0.960289856497536*[1. -1.]];
        weights(ii) = [0.362683783378362*[1.  1.] 0.313706645877887*[1.  1.] 0.222381034453374*[1.  1.] 0.101228536290376*[1.  1.]];
    case (9)
        ii          = [5,6,4,7,3,8,2,9,1];
        nodes(ii)   = [0.                         0.324253423403809*[1. -1.] 0.613371432700590*[1. -1.] 0.836031107326636*[1. -1.] 0.968160239507626*[1. -1.]];
        weights(ii) = [0.330239355001260          0.312347077040003*[1.  1.] 0.260610696402935*[1.  1.] 0.180648160694857*[1.  1.] 0.081274388361574*[1.  1.]];
    case (10)
        ii          = [6,5,7,4,8,3,9,2,10,1];
        nodes(ii)   = [0.148874338981631*[1. -1.] 0.433395394129247*[1. -1.] 0.679409568299024*[1. -1.] 0.865063366688985*[1. -1.] 0.973906528517172*[1. -1.]];
        weights(ii) = [0.295524224714753*[1.  1.] 0.269266719309996*[1.  1.] 0.219086362515982*[1.  1.] 0.149451349150581*[1.  1.] 0.066671344308688*[1.  1.]];
    otherwise
        error('Gaussian quadrature is not implemented for this number of points: %i.', NPoints)
end

numEle      = length(Domain)-1;
%[PTS,WTS]   = gaussianQuadratureData([-1,1],'Degree',Degree);
%numPts      = length(WTS);
Pts  = zeros(1,NPoints*numEle);
Wts = zeros(1,NPoints*numEle);

for j = 1:numEle % Loop over elements
    ind                    = (j-1)*NPoints + (1:NPoints);
    dom                    = Domain(j+(0:1));
    Pts(ind)  = 0.5*( diff(dom)*nodes + sum(dom) );
    Wts(ind) = 0.5*diff(dom).*weights;
end

% Scale and translate function arguments, and scale weights for this domain.
% Pts = 0.5*( diff(Domain)*nodes + sum(Domain) );
% Wts = 0.5*diff(Domain).*weights;

end
