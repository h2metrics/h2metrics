function res = pVectors(pts,tri)
% computes (a representation of) the p-vector (ie tangent vector for curve and normals for surfaces). 
%
% Input:
%  pts: list of vertices (n x d matrix)
%  tri: list of edges (T x M matrix of indices)
%
% Output:
%  res: list of p-vectors (d x M matrix)
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)

M = size(tri,2);
d = size(pts,2);

if (M==2)
    
    res=pts(tri(:,2),:)-pts(tri(:,1),:);
    
elseif (M==3) && (d==3) % surfaces

    u=pts(tri(:,2),:)-pts(tri(:,1),:); 
    v=pts(tri(:,3),:)-pts(tri(:,1),:);
    
    res =[u(:,2).*v(:,3)-u(:,3).*v(:,2),...
          u(:,3).*v(:,1)-u(:,1).*v(:,3),...
          u(:,1).*v(:,2)-u(:,2).*v(:,1)];
end

end