%% createTensorCollocationMatrix
%
% Computes collocation matrix to evaluate a path or its derivatives
% on a cartesian product of evaluation sites
%
% Input
%   evalSitesS, evalSitesT
%       Positions, where to evaluate the path
%   derS, derT
%       Which derivatives to evaluate; uses the same convention as spcol
%       To evaluate function use derS=1, to evaluate first derivative
%       use derS=2, etc.
%   splineData
%       General information about the splines used.
%
% Output
%   B
%       The collocation matrix
%
function B = createTensorCollocationMatrix( evalSitesS, evalSitesT, ...
    derS, derT, splineData )

N = splineData.N;
nS = splineData.nS;
Nt = splineData.Nt; 
nT = splineData.nT;

noEvalSitesS = length(evalSitesS);
noEvalSitesT = length(evalSitesT);

knotsS = splineData.knotsS;
knotsT = splineData.knotsT;
    
B_S = spcol( knotsS, nS+1, ...
              brk2knt( evalSitesS, derS ), 'sparse');
B_S = [ B_S(:,1:nS) + B_S(:,end-nS+1:end), ... % For periodic spline
        B_S(:,nS+1:end-nS) ];
B_S = B_S(derS:derS:end,:); % Take only the derivative we want

B_T = spcol( knotsT, nT+1, ...
             brk2knt( evalSitesT, derT ), 'sparse');              
B_T = B_T(derT:derT:end,:);

nnzmaxST = noEvalSitesS*(nS+1) * noEvalSitesT*(nT+1);
B = spalloc(noEvalSitesS*noEvalSitesT, N*Nt, nnzmaxST); % Preallocate

for jj = 1:Nt
    for ii = 1:N
        B(:,ii+(jj-1)*N) = ...
            reshape(B_S(:,ii)*B_T(:,jj)', noEvalSitesS*noEvalSitesT, 1);
    end
end

end