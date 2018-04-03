function easyfun_driver(n1,n2)
% a simple test of HANSO
fprintf('First easy test problem: HANSO calls BFGS starting from 10 randomly generated starting points\n\n')
pars = pars_easyfun(n1,n2);
[x,f] = hanso(pars); % or hanso(pars,options): see documentation (type "help hanso")
%
fprintf('\nSecond easy test problem: HANSO calls BFGS starting from 10 randomly generated starting points\n\n')
pars = pars_easyfun2(n1,n2);
[x,f] = hanso(pars); % or hanso(pars,options): see documentation (type "help hanso")