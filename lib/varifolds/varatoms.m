function [X,Xi] = varatoms(pts,tri)
% [X,tf,Xi] = FCATOMS(pts,f,tri,signal_type) compute the dirac representation of meshes (1D or 2D)
%
% Input :
%  pts :  matrix with vertices (one column per dimension, matrix nxd)
%  tri : connectivity matrix. if tri is of size (M x 2) this is a 1d current 
%  	(so each line of tri contains the two indexes of the points defining 
%  	an oriented segment), if tri is (M x 3) this to be considered 
%  	as a 2d current (so each line of tri contains the three indexes of 
%  	the points defining an oriented triangle).
%  f : column vector of functional values (n x1)
%
% Output
% X : the matrix of the centers of the faces   (M x d)
% Xi: is the matric of p-vectors (tangent (1d current) or normal 2d current) (M x 2 or 3)
%
% To do : generalize to measure and simplexes
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, G. Nardi, A. Trouve (2012-2015)


[T,M]=size(tri);
d = size(pts,2);

%normals 
Xi = pVectors(pts,tri) / (M-1); % normals or tangent vectors

%centers
X=sum(reshape(pts(tri(:),:)',d,T,M ),3)'/M; % center of the faces

end

