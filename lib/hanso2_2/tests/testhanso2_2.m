function testhanso2_2(wantpause)
if nargin < 1
    wantpause = 0;
end
fprintf(' test suite for HANSO 2.2\n')
fprintf(' mostly same as suite for HANSO 2.0 except that we now include the easyfun and easyfun2 examples\n')
fprintf(' AND we test limited memory BFGS by invoking it explicitly at end\n')
fprintf(' AND we test gradient sampling by invoking it explicitly at end\n')
fprintf(' because the new defaults in HANSO 2.1 and 2.2 are\n')
fprintf(' respectively NO limited memory BFGS and NO gradient sampling\n\n')
fprintf(' default is 10 BFGS starting points, options.normtol 1e-4, options.maxit 1000\n')
fprintf(' but we also include other choices for normtol, maxit in some of the tests below\n\n ')


% start with easyfun and easyfun2, which were previously not included in the test suite for HANSO 2.1
fprintf(' First easy test problem\n')
pars = pars_easyfun(10,10);
[x,f] = hanso(pars); % or hanso(pars,options): see documentation (type "help hanso")

if wantpause, fprintf('\n\nhit any key to continue\n'), pause, end
fprintf(' Second easy test problem\n')
pars = pars_easyfun2(10,10);
[x,f] = hanso(pars); 

if wantpause, fprintf('\n\nhit any key to continue\n'), pause, end
fprintf(' Nesterov SMOOTH Chebyshev-Rosenbrock function\n')
fprintf(' this is very hard so sometimes the optimality termination is satisfied\n')
fprintf(' not by norm of gradient being small but norm of smallest vector in convex\n')
fprintf(' hull of gradients small, even though f is actually smooth\n')
pars = pars_yurirosen(7)
[x, f, loc, X, G, w, H] = hanso(pars); % default options

if wantpause, fprintf('\n\nhit any key to continue\n'), pause, end
fprintf(' reduce tolerance and increase maximum iterations for BFGS\n')
options.normtol = 1e-6;
options.maxit = 10000;
[x, f, loc, X, G, w, H] = hanso(pars,options);
clear options

if wantpause, fprintf('\n\nhit any key to continue\n'), pause, end
fprintf(' Nesterov FIRST NONSMOOTH Chebyshev-Rosenbrock function\n')
fprintf(' this is partly smooth and the only Clarke stationary point is the local\n')
fprintf(' minimizer, but the problem is tough even for n=3\n')
pars = pars_yurirosen_ns1(3)
[x, f, loc, X, G, w, H] = hanso(pars); % default options

if wantpause, fprintf('\n\nhit any key to continue\n'), pause, end
fprintf(' Nesterov SECOND NONSMOOTH Chebyshev-Rosenbrock function\n')
fprintf(' this is not regular and hence has multiple Clarke stationary points, note\n')
fprintf(' how different starting points result in termination at any one of 4\n')
fprintf(' stationary points for n=3\n')
pars = pars_yurirosen_ns2(3)
[x, f, loc, X, G, w, H] = hanso(pars); % default options

if wantpause, fprintf('\n\nhit any key to continue\n'), pause, end
fprintf(' Anstreicher-Lee eigenvalue product - 25 variables (nonsmooth)\n')
% in the call to parsdef, the first 5 refers to the size of the matrix and
% the second 5 to the rank: this is better than an alternative way of posing
% the full rank problem as explained in the comments
pars = parsdef(5,5,1) % 1 (jon36), not 0, as 5,5,0 is a smooth problem
[x, f, loc, X, G, w, H] = hanso(pars); % default options

if wantpause, fprintf('\n\nhit any key to continue\n'), pause, end
fprintf(' Anstreicher-Lee eigenvalue product - 25 variables again, but this time invoke gradient sampling\n')
fprintf(' which is slower, but has convergence guarantees\n')
pars = parsdef(5,5,1);
% give BFGS small maxit and demanding tolerance so Gradient Sampling has a chance to improve it
options.maxit = 50;
options.normtol = 1e-6;
options.samprad = [1e-3 1e-4 1e-5]; % as default for options.evaldist is 1e-4
[x, f, loc, X, G, w, H] = hanso(pars,options);
clear options

if wantpause, fprintf('\n\nhit any key to continue\n'), pause, end
fprintf(' Anstreicher-Lee eigenvalue product - 100 variables (nonsmooth)\n')
pars = parsdef(10,10,0)  % 0 is the default: kurt63
[x, f, loc, X, G, w, H] = hanso(pars); % default options

if wantpause, fprintf('\n\nhit any key to continue\n'), pause, end
fprintf(' Anstreicher-Lee eigenvalue product - 400 variables (nonsmooth)')
fprintf(' takes longer, so use only one starting point\n')
pars = parsdef(20,20,0) % 0 is the default: kurt63
options.x0 = randn(pars.nvar,1);
[x, f, loc, X, G, w, H] = hanso(pars,options);
% DO NOT clear options

if wantpause, fprintf('\n\nhit any key to continue\n'), pause, end
fprintf(' Anstreicher-Lee eigenvalue product - 400 variables again, but this time invoke limited memory BFGS\n')
fprintf(' which is faster, but does not get such a good answer\n')
fprintf(' use same x0 as full BFGS run') 
pars = parsdef(20,20,0) % 0 is the default: kurt63
options.nvec = 3; % limited memory size: maintain 5 BFGS updates
[x, f, loc, X, G, w, H] = hanso(pars,options);
clear options



