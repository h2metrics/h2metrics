function [pars,options] = pars_yurirosen_ns1(n,scalex0)
% Nesterov's Chebyshev-Rosenbrock nonsmooth version 1
pars.fgname = 'yurirosen_ns1';
pars.title = 'Nesterov-Chebyshev-Rosenbrock Nonsmooth#1'; 
pars.nvar = n;
pars.optval = 0;
pars.title = 'Nesterov-Chebyshev-Rosenbrock 1';
if nargin >= 2 % if not, don't set options.x0
    options.x0 = scalex0*ones(n,1);  % BFGS fails immediately if scalex0 is 1, but GS does not
    options.x0(1) = -scalex0;
end
options.H0 = eye(n); % so does not do initial scaling for BFGS
options.fvalquit = 1e-15;
options.maxit = 10^(n+1); % !!!  for BFGS
% options.phasenum = [1 10 6]; % for HANSO
% options.phasemaxit = [10^n 10^n 1000]; % for HANSO
