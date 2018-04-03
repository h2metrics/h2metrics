function pars = parsdef(n, r, Adata, K, rho)
% call:  pars = parsdef(n, r, Adata, K, rho)
% set up an Anstreicher-Lee log of eigenvalue products problem.
% n is dimension of A 
% r = 0 means full rank problem, implemented with X the variable, diag(X)=1,
% otherwise low rank with prescribed rank r, that is X=VV', V is n by r
% interestingly, turns out to be best to use r=n, in preference to r=0
% Adata=0 means Kurt's original 63-dim data, Adata=1 is Jon's subsequent
% 36-dim data, Adata=-1 means randomly generated, positive definite A
% K is the number of eigenvalues in the product: default [n/2]
% rho is the penalty parameter, only makes sense if r=0

% parsdef replaces parsdef_lowrank and parsdef_fullrank.
if nargin < 5
    rho = 1;
end
if nargin < 4
    K = floor(n/2);
end
if nargin < 3
    Adata = 0;  % Kurt 63: the one we used for Gradient Sampling paper
end
if round(r) ~= r
    error('r must be an integer')
end
if Adata == 1
   load jon36   %already symmetric 36 by 36
   A = jon36;
   pars.data = 'jon36';
elseif Adata == 0
   load kurt63
   A = reshape(kurt63',63,63);  % note the transpose
   pars.data = 'kurt63';
else
   B = randn(n);
   A = B'*B;
   A = 0.5*(A+A');
   pars.data = 'rand';
end
A = A(1:n,1:n); 
A = A/max(max(abs(A))); % so largest entry is 1
pars.A = A;
pars.n = n;
pars.K = K; 
% pars.nonneg = 0;  % no nonnegativity constraint: was in obsolete code eigprod_gen_slow.m
if r == 0 % full rank problem
    pars.formulation = 1;
    pars.nvar = n*(n-1)/2;
    pars.indxlower = trilindx(n);  % to work with the new buildmat.m
    pars.rho = rho;  % impose pos def constraint with exact penalty
                     % even with modest pars.rho, this seems to create
                     % trouble with nonsmoothness that is not present in
                     % the eigenvalue product
    pars.fgname = 'logeigprod'; % LOG of eigenvalue product
    pars.title = sprintf('Log Eig Prod, Data=%s, N=%d, full rank, K=%d, nvar=%d', pars.data, n, pars.K, n*(n-1)/2);
else
    pars.nvar = n*r;
    pars.rank = r;
    pars.fgname = 'logeigprodlowrank'; % LOG of eigenvalue product, low rank formulation
    if Adata == 1 % distinguish the cases without a weird title
        pars.title = sprintf('LOG EIG PROD, N=%d, r=%d, K=%d, n=%d', n, r, pars.K, n*r);
        pars.varytitle = 'LOG EIG PRODUCT';
    else
        pars.title = sprintf('Log eig prod, N=%d, r=%d, K=%d, n=%d', n, r, pars.K, n*r);
        pars.varytitle = 'log eig product';
    end
end