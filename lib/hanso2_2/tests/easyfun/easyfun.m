function [f,g]=easyfun(x,pars)
% an easy example for testing HANSO or BFGS
% the first n1 terms of the objective are nonsmooth at 0 and the last n2 are smooth (quadratic)
% f = sum from i=1 to n1 of abs(x(i)) + sum from i=n1+1 to n1+n2 of 0.5*x(i)^2
% this example is convex so it can be solved more efficiently by other methods
% the solution is x=[0,...,0]'
n1 = pars.n1;
n2 = pars.n2;
f = sum(abs(x(1:n1))) + 0.5*sum(x(n1 + 1: n1 + n2).^2);
g1 = sign(x(1:n1)); 
g2 = x(n1 + 1 : n1 + n2);
g = [g1; g2];