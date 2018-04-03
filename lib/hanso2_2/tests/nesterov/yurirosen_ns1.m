function [f,g] = yurirosen_ns1(x,pars)
n = pars.nvar;
g = zeros(n,1);
% this is version 1 of Yuri's nonsmooth Chebyshev Rosenbrock, but there are
% two versions of this, according to the choice of the initial term
% f = abs(1-x(1))/4;       % in this case the dimension of the V space is 2
% g(1) = sign(x(1)-1)/4;   % (not sure if f is regular at the solution)
                   % BFGS does well for n=4 50% of time, and occasionally for n=5
% we go with the 2nd variant of version 1 in the paper, because it
% illustrates the difficulty BFGS can have when there is a nontrivial
% U space and it takes a long time to traverse along it
f = (1-x(1))^2/4;  % in this case the U and V spaces both have dim 1
g(1) = -(1-x(1))/2; % problem is much more difficult: cannot solve n=4
                    % accurately, probabaly because of rounding
    
for i=1:n-1
    f = f + abs(1+x(i+1)-2*x(i)^2);
    r = sign(1+x(i+1)-2*x(i)^2);
    g(i+1) = g(i+1) + r;
    g(i) = g(i) - 4*x(i)*r;
end