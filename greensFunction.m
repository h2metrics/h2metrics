%% greensFunction
% Computes the Green's function associated to a metric
%
% NOTE: First version, signature might change.
%
% NOTE: Only works if all roots of characteristic polynomial are distinct
%
% Input
%   a
%       Coefficient vector [L2 H1 H2 ... Hn]
%   ell
%       Length of the curve; assumed to be constant speed
%   evalPts
%       Where to evaluate the function; should be in [0, 2*pi]
%   y0
%       Where to center the Green's function; G(. - y0)
%
% Output
%   y
%       Green's function evaluated at evalPts; y has same shape as evalPts
function y = greensFunction(a, ell, evalPts, y0)

n = length(a) - 1;
b = zeros([2*n+1,1]);
ga = 2*pi/ell;
for jj = 1:n+1
    b(2*jj-1) = a(jj) * ga;
    ga = ga * (-1)*(ell/(2*pi))^2;
end

y = greensConstCoeff(b, evalPts, y0);

% n = length(a) - 1;
% 
% % Normalize coefficients and setup polynomial
% a2 = zeros([n+1, 1]);
% for jj = 1:n+1
%     a2(jj) = a(jj) * (-1)^(jj-1) * (ell/(2*pi))^(-2*jj+3);
% end
% a2 = a2 / a2(n+1);
% 
% p = zeros([1, 2*n+1]);
% for jj = 1:n+1
%     p(2*jj-1) = a2(jj);
% end
% p = flip(p); % The polynomial
% 
% la = roots(p);
% 
% A = zeros([2*n, 2*n]);
% for jj = 1:2*n
%     A(jj,:) = la.'.^(jj-1); % Note the non-conjugating transpose here
% end
% 
% b = zeros([2*n, 1]);
% b(2*n) = 1;
% 
% c = A \ b;
% 
% be = exp(y0*la) - exp((2*pi+y0)*la);
% b = c ./ be;
% 
% evalPts = mod(evalPts-y0, 2*pi) + y0;
% 
% y = zeros(size(evalPts));
% for jj = 1:(2*n)
%     y = y + b(jj) * exp(la(jj)*evalPts);
% end
% y = y / a2(n+1);
% y = real(y);

function y = greensConstCoeff(b, evalPts, y0)

n = length(b) - 1;

p = zeros([1, n+1]);
for jj = 1:n+1
    p(jj) = b(jj);
end
p = flip(p); % The polynomial

la = roots(p); % We assume all roots are distinct

A = zeros([n, n]);
for jj = 1:n
    A(jj,:) = la.'.^(jj-1); % Note the non-conjugating transpose here
end

d = zeros([n, 1]); % rhs
d(n) = 1;

c = A \ d;

be = exp(y0*la) - exp((2*pi+y0)*la);
c = c ./ be;

evalPts = mod(evalPts-y0, 2*pi) + y0; % lying in [y0, y0+2*pi]

y = zeros(size(evalPts));
for jj = 1:n
    y = y + c(jj) * exp(la(jj)*evalPts);
end
y = real(y);

