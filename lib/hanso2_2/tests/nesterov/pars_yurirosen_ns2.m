function [pars,options] = pars_yurirosen_ns2(n,scalex0)
% Nesterov's Chebyshev-Rosenbrock nonsmooth version 2
pars.fgname = 'yurirosen_ns2';
pars.title= 'Nesterov-Chebyshev-Rosenbrock Nonsmooth#2';
%%%% NOTE THE 2, MORE INTERESING OF THE TWO,
         % because for this version, there are 2^(n-1) Clarke stationary
         % points, all of which are attractors for BFGS
pars.nvar = n;
pars.optval = 0;
if nargin >= 2 % otherwise don't set options.x0
    options.x0 = scalex0*ones(n,1);  % BFGS fails immediately if scalex0 is 1, but GS does not
    options.x0(1) = -scalex0;
end
options.H0 = eye(n); % so does not do initial scaling for BFGS
options.maxit = 10^(n+1); % !!!  for BFGS
options.fvalquit = 1e-15;
% options.phasenum = [1 10 6]; % for HANSO
% options.phasemaxit = [10^n 10^n 1000]; % for HANSO
pars.title = 'Nesterov-Chebyshev-Rosenbrock'; % omit the 2
pars.varytitle = pars.title;
