function [f,g,B] = logeigprodlowrank(x,pars)
% objective is LOG of product of largest pars.K eigenvalues as coded
% in eigprodlowrank.  For large n, the product gets too small so 
% we need to take logs.
n_hardwire = 63;
r_hardwire = 5;
data_hardwire = 'kurt63';
if nargin < 2 % when using PBUN, nargin can only be 1
    n = n_hardwire;  % hardwire size of problem and the data
    r = r_hardwire;
    K = floor(n_hardwire/2);
    if ~strcmp(data_hardwire,'kurt63')
        load jon36   %already symmetric 36 by 36
        A = jon36;
    else
        load kurt63;
        A = reshape(kurt63',63,63);  % note the transpose
    end
    A = A(1:n,1:n); 
    A = A/max(max(abs(A))); % so largest entry is 1
    pars.A = A;
    pars.n = n;
    pars.K = K; 
    pars.nvar = n*r;
    pars.rank = r;
else % uncomment this if comparing our codes with PBUN
%     if ~strcmp(pars.data,data_hardwire) | pars.n ~= n_hardwire | pars.rank ~= r_hardwire
%         fprintf('pars parameters do not match hardwired values\n')
%         keyboard
%     end
end
[fp, gp, B] = eigprodlowrank(x,pars); % B is the matrix
f = log(fp);
g = gp/fp;