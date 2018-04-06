%% tangentPCA
%
% Visualizes PCA in two dimensions
%   U
%       Basis of Eigenvectors
%   lambda 
%       Eigenvalues
%   v
%       List of initial velocities of the geodesics from the 
%       reference curve to the other curves
%
% Output
% Plot of a 2-d projection of the data.
function visualizePCA(U,Lambda,G,vList,splineData)
    noCurves = length(vList);
    N=splineData.N;
    V = U(:,1:2);
    rootG = sqrtm(G);
    pts2d = zeros([2, noCurves]);
    for jj = noCurves:-1:1
        v = reshape(vList{jj}, [N*2, 1]);
        pts2d(:,jj) = V' * rootG * v;
        pts2d(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
    end
    for i=1:noCurves
    plot(pts2d(1,i),pts2d(2,i),'bx')
    hold on
    end 