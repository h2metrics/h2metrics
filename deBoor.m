%% deBoor
%
% Uses de Boor algorithm to evaluate the spline with given knots and
% control points d at given points.
%
% IMPORTANT
%   For non-periodic input the function is only correct up to and
%   including the first derivative. For higher derivatives boundary knots
%   create problems.
%
% Input
%   knots
%       The knot vector.
%       Only uniformly spaced inner knots are supported.
%   nS
%       Order of the spline. Note that this corresponds to nS+1 for spcol.
%   d
%       Control points of the spline. If d is a matrix, we evaluate the
%       splines corresponding to each column of d.
%   evalPoints
%       Points, where to evaluate the spline. Has to be a column vector.
%   order
%       How many derivatives do we want to evaluate. To evaluate the first
%       derivative, set
%           order=2
%       Only up to three derivatives (order={1,2,3,4}) are supported 
%       at the moment.
%
% Optional parameters
%   'periodic' = {true, false (default)}
%       If true, then the vector d contains the independent control points
%       of a periodic spline.
%
% Output
%   y
%       Column vector or matrix containing the spline and its derivatives
%       at the given points. Dimensions of y are
%           [length(evalPoints)*order, size(d,2)]
%       The function evaluations are located at
%           y(1:order:end,:)
%       and first derivatives at
%           y(2:order:end,:)
%       This corresponds to the output of spcol.
function y = deBoor( knots, nS, d, evalPoints, order, varargin )

periodic = false;

ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'periodic'
                ii = ii + 1;
                if isa(varargin{ii}, 'integer') || ...
                        isa(varargin{ii}, 'logical')
                    periodic = logical(varargin{ii});
                else
                    error('Invalid value for option ''periodic''.');
                end
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1; 
    end
end

if order > 4
    error ('Order too high. Only three derivatives supported.');
end

noPoints = length(evalPoints);
noCols = size(d, 2);
N = length(knots) - 2*nS - 1;
y = zeros([noPoints*order, noCols]);

t_min = knots(nS+1); % Limits for inner knots
t_max = knots(end-nS);
t_diff = t_max - t_min;

t_work = evalPoints;
if size(evalPoints,1) == 1 % evalPoints has to be a column vector
    t_work = evalPoints';
end

if periodic
    d_nonper = [ d; d(1:nS,:) ];
else
    d_nonper = d;
end

%% de Boor's algorithm vectorized
r = floor((t_work-t_min)*N / ...
        (t_max-t_min)+nS)+1; % Here we assume uniform spacing of knots.
r = r - (r == N + nS + 1); % Takes care of the case t==t_max.
s = zeros(noPoints, noCols, 2*nS);
for jj = 1:noCols
    for kk = 1:2*nS
        s(:,jj,kk) = knots(r-nS+kk);
    end
end

d_work = zeros(noPoints,noCols,nS+1);
for jj=1:noCols
    for kk = 1:nS+1
        d_work(:,jj,kk) = d_nonper(r-nS-1+kk, jj);
    end
end

al = zeros(noPoints, noCols, nS);
for kk = 1:nS
    % al(:,kk:nS) = (evalPoints*ones([1, nS-kk+1]) - ...
    %   s(:,kk:nS)) ./ (s(:,nS+1:2*nS-kk+1) - s(:,kk:nS));
    for jj = kk:nS
        al(:,:,jj) = (t_work*ones([1, noCols]) - s(:,:,jj)) ./ ...
            (s(:,:,nS+jj-kk+1) - s(:,:,jj));
    end
    
    % d_work(:,:,kk+1:nS+1) = (1-al(:,:,kk:nS)) .* d_work(:,:,kk:nS) + ...
    %   al(:,:,kk:nS) .* d_work(:,:,kk+1:nS+1);
    for jj = nS+1:-1:kk+1
        d_work(:,:,jj) = (1-al(:,:,jj-1)) .* d_work(:,:,jj-1) + ...
            al(:,:,jj-1) .* d_work(:,:,jj);
    end
    
    if order >= 4 && kk == nS-3
        y(4:order:end,:) = 1/6 * nS*(nS-1)*(nS-2) * ...
            (d_work(:,:,kk+4) - 3*d_work(:,:,kk+3) + ...
             3*d_work(:,:,kk+2) - d_work(:,:,kk+1)) ./ (t_diff/N)^3;
    end
    
    if order >= 3 && kk == nS-2
        y(3:order:end,:) = 0.5 * nS*(nS-1) * ...
            (d_work(:,:,kk+3) - 2*d_work(:,:,kk+2) + ...
             d_work(:,:,kk+1)) ./ (t_diff/N)^2;
    end
    
    if order >= 2 && kk == nS-1
        y(2:order:end,:) = nS * (d_work(:,:,kk+2) - d_work(:,:,kk+1)) ...
            ./ (t_diff/N); % ./ (s(:,:,nS+1) - s(:,:,nS));
    end
end
y(1:order:end,:) = d_work(:,:,nS+1);

end