function [P, Pinv, PBlock, PinvBlock] = ...
    createPreconditioner(lambda, splineData)

%% Create single matrix
b1 = sqrt(splineData.a(1) * lambda);
b2 = sqrt(splineData.a(3) / lambda^3);

N = splineData.N;
dSpace = splineData.dSpace;

quadData = splineData.quadData;
noQuadSites = quadData.noQuadPointsS;
quadWeights = quadData.quadWeightsS;

B = quadData.B_S;
Buu = quadData.Buu_S;

sparseW = sparse(1:noQuadSites, 1:noQuadSites, quadWeights);

A1 = B' * sparseW * B;
A2 = B' * sparseW * (b1*B - b2*Buu);

% P = A1 \ A2;
Pinv = A2 \ A1;
Pinv = full(Pinv);

% Create some sparsity by setting consecutive diagonals equal to 0
% symmetry: (1,2) - (n,1), (1,3) - (n-1,1), ...
% k = n/2-1
% diagonal is (1,1+k+1), (2,1+k+j), ..., (jm, 1+k+jm)
n = size(Pinv, 1);
kk1 = floor(n/2) - 1;
kk2 = ceil(n/2) - 1;
while kk1 >= 1 && kk2 <= n && ...
      nnz(Pinv) > numel(Pinv) * splineData.precondSparsity
    
    for jj = 1:n-kk1-1
      Pinv(jj, 1+kk1+jj) = 0;
      Pinv(1+kk1+jj, jj) = 0;
    end
    
    for jj = 1:n-kk2-1
      Pinv(jj, 1+kk2+jj) = 0;
      Pinv(1+kk2+jj, jj) = 0;
    end
  
    kk1 = kk1 - 1;
    kk2 = kk2 + 1;
end

Pinv = sparse(Pinv);
P = inv(Pinv);

%% Create block matrix for path
Nt = splineData.Nt;

[row, col, val] = find(P);
newRow = repmat(row,1,(Nt-2)*dSpace) + ...
            repmat(0:(Nt-2)*dSpace-1,length(row),1)*N;
newCol = repmat(col,1,(Nt-2)*dSpace) + ...
            repmat(0:(Nt-2)*dSpace-1,length(col),1)*N;
newVal = repmat(val,1,(Nt-2)*dSpace);

PBlock = sparse(newRow, newCol, newVal);


[row, col, val] = find(Pinv);
newRow = repmat(row,1,(Nt-2)*dSpace) + ...
            repmat(0:(Nt-2)*dSpace-1,length(row),1)*N;
newCol = repmat(col,1,(Nt-2)*dSpace) + ...
            repmat(0:(Nt-2)*dSpace-1,length(col),1)*N;
newVal = repmat(val,1,(Nt-2)*dSpace);

PinvBlock = sparse(newRow, newCol, newVal);


