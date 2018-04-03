function [f,g,B,Vhat] = eigprodlowrank(v, pars)
% function [f,g,B,Vhat,DVhat,W,D] = eigprodlowrank(v, pars)
% a problem from Kurt Anstreicher and Jon Lee: low rank version.
% Reformulated Feb and Apr 2005.  New code is much faster than
% either of the original formulations. See notes of Feb 14 at ibm.
% Most important change was avoiding inner loop in gradient calculation.
%
% objective is product of largest pars.K eigenvalues of A o Xhat where A 
% (passed as pars.A) is fixed symmetric matrix, X = Vhat Vhat', where V 
% is pars.n by pars.rank, whose columns are passed in the vector v,
% and Vhat is obtained from V by Vhat(i,:)=V(i,:)/norm(V(i,:).
% By construction Xhat is symmetric, positive semidefinite, low-rank, 
% and has diagonal all ones, so there is no need for a penalty term.
% 
% Returns function f, gradient g and also the matrix B whose eigenvalue
% product is being computed
if isfield(pars, 'formulation')
    disp('there is no longer a formulation field in low rank code')
end
m = pars.nvar;  % expensive to evaluate these fields inside loops
n = pars.n;
A = pars.A;  
if size(A,1) ~= n | size(A,2) ~= n
   error('mismatch between pars.A and pars.n')
end
K = pars.K;
r = pars.rank;
% construct V from v and Vhat from V
V = reshape(v, n, r);
for i=1:n
    v = V(i,:);
    Vhat(i,:) = v/norm(v);
end
W=Vhat*Vhat';
B = A .* (Vhat*Vhat');   % Hadamard product A o Xhat
[Q, evals] = eig(B);   
evals = diag(evals);           % vector of eigenvalues
[junk, indx] = sort(-evals);   % sort descending
lambda = evals(indx);          % sorted eigenvalues
Q = Q(:,indx);                 % and corresponding eigenvectors
eprod = prod(lambda(1:K));         % objective function
f = eprod;
% now the gradient
I = eye(n,r);
Ibig = eye(n);
% preconstruct the lambda products
for p = 1:K
   lambdaprod(p) = eprod/lambda(p);
end
Q = Q(:,1:K);
QlambdaprodQt = Q*diag(lambdaprod)*Q';
for i = 1:n     
   for j = 1:r                 % gradient with respect to V(i,j)
      ej = I(j,:);
      v = V(i,:);
      nmv = norm(v);
      w = ((ej - v*(V(i,j)/nmv^2)))'/nmv;
      % eibig = Ibig(:,i);
      % DVhat = eibig*w';
      % D = eibig*w'*Vhat';
      % D = D + D';
      % M = A .* D;
      what = Vhat*w; % pronounced w-hat!
      z = A(:,i) .* what;
      s = 0;  
      % the following inner loop is very slow
      % for p = 1:K              % run through the eigenvalues in the product
         % this version is slowest of all
         % s = s + (eprod/lambda(p))*Q(:,p)'*M*Q(:,p);
         % this is not as bad, but still not good
         % q = Q(:,p);
         % s = s + 2*(eprod/lambda(p))*q(i)*(q'*z);
      %end
      eibig = Ibig(:,i);
      s = 2*eibig'*QlambdaprodQt*z;
      G(i,j) = s;
  end
end
g = G(:);   % turn from matrix into vector
