function [f,g] = yurirosen_ns2(x,pars)
 %%%% NOTE THE 2, MORE INTERESING OF THE TWO,
 % because for this version, there are 2^(n-1) Clarke stationary
 % points, all of which are attractors for BFGS
n = pars.nvar;
f = abs(1-x(1))/4;
g = zeros(n,1);
g(1) = sign(x(1)-1)/4;
for i=1:n-1
    f = f + abs(1+x(i+1)-2*abs(x(i)));
    r = sign(1+x(i+1)-2*abs(x(i)));
    g(i+1) = g(i+1) + r;
    g(i) = g(i) - 2*sign(x(i))*r;
end