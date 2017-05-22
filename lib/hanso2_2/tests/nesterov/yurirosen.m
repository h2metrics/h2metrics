function [f,g] = yurirosen(x,pars)
% based on Chebyshev polynomials
% x_{i+1} = 2x_i^2 - 1 = T_2(x_i) = T_{2^i}(x_1) = cos(2^i arccos(x_1))
n = pars.nvar;
f = (1-x(1))^2/4; % the 1/4 gives initial value 1 and prevents method
%                   skipping over the nasty prescribed path
g = zeros(n,1);
g(1) = (x(1)-1)/2;
for i=1:n-1
    f = f + (1+x(i+1)-2*x(i)^2)^2;
    r = 1+x(i+1)-2*x(i)^2;
    g(i+1) = g(i+1) + 2*r;
    g(i) = g(i) - 8*x(i)*r;
end