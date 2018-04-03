function [f,g]=easyfun2(x,pars)
% a second easy example for testing HANSO or BFGS
% the first part of the objective is the 2-norm of the first n1 variables
% which is nonsmooth at 0, and the last part is smooth (the norm squared)
% f = sqrt(sum from i=1 to n1 of x(i)^2) + sum from i=n1+1 to n1+n2 of x(i)^2
% this example is convex so it can be solved more efficiently by other methods
% the solution is x=[0,...,0]'
n1 = pars.n1;
n2 = pars.n2;
f =  sqrt(sum(x(1:n1).^2)) + sum(x(n1 + 1: n1 + n2).^2);
g1 = x(1:n1)/sqrt(sum(x(1:n1).^2));
g2 = x(n1 + 1 : n1 + n2);
g = [g1; g2];