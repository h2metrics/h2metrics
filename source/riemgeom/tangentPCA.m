%% tangentPCA
%
% Principal component analysis of a set of curves at the tangent space of a
% reference curve
% Input
%   dList
%       List of curves
%   d0  reference curve
%   splineData
%       General information about the splines used.
%
% Output
%   U
%       Basis of Eigenvectors
%   lambda 
%       Eigenvalues
%   v
%       List of initial velocities of the geodesics from the 
%       reference curve to the other curves
%   G
%       metric matrix at d0
function [U,Lambda,G,vList] = TangentPCA(dList,d0,splineData)
    noCurves = length(dList);
    G = metricMatrixH2(d0, splineData);
    rootG = sqrtm(G);
    vList = {};
    for jj = noCurves:-1:1
        [~, optPath, ~, ~]= geodesicBvp(d0, dList{jj},splineData ,splineData.options);
        vList{jj} = pathVelocity(optPath, 0, splineData);
    end
    N=splineData.N;
    Sigma = zeros([N*2, N*2]);
    for jj = noCurves:-1:1
        v = reshape(vList{jj}, [N*2, 1]);
        Sigma = Sigma + v * v';
    end
    Sigma = 1./(noCurves-1) * rootG * Sigma * rootG;
    [U, Lambda] = eig(Sigma);
    Lambda = real(diag(Lambda));
    
