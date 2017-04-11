function [dxg2]=chain_rule(x,tri,DXg,DXig)
%
% Computation of the gradient with respect to x by distributing the previous gradients
% on points of the initial shape. (Chain's rule). 
%
% Input:
%  x: list of vertices (n x d matrix)
%  G: list of edges (T x M matrix)
%  DXg: derivative wrt X (center of the cells: T x d matrix)
%  DXig: derivative wrt Xi (p-vectors: T x d matrix)
%
% Output:
%  dxg: derivative wrt x (n x d matrix)
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)


[nx,d]=size(x);
[~,M] = size(tri);


if d==2
    U2 =  [accumarray(tri(:),repmat(DXg(:,1),M,1),[nx,1],[],0),...
           accumarray(tri(:),repmat(DXg(:,2),M,1),[nx,1],[],0)] / M;
elseif d==3
    U2 =  [accumarray(tri(:),repmat(DXg(:,1),M,1),[nx,1],[],0),...
           accumarray(tri(:),repmat(DXg(:,2),M,1),[nx,1],[],0),...
           accumarray(tri(:),repmat(DXg(:,3),M,1),[nx,1],[],0)] / M;
end


if M==2 % curve case
    if (d==2)
        dxg2 = U2 + [accumarray(tri(:),[-DXig(:,1);DXig(:,1)] ,[nx,1],[],0),...
                     accumarray(tri(:),[-DXig(:,2);DXig(:,2)] ,[nx,1],[],0) ];
    elseif (d==3)
        dxg2 = U2 + [accumarray(tri(:),[-DXig(:,1);DXig(:,1)] ,[nx,1],[],0),...
                     accumarray(tri(:),[-DXig(:,2);DXig(:,2)] ,[nx,1],[],0),...
                     accumarray(tri(:),[-DXig(:,3);DXig(:,3)] ,[nx,1],[],0) ];
    end
    
elseif (M==3) && (d==3) % surface case
    
    Xa=x(tri(:,1),:);
    Xb=x(tri(:,2),:);
    Xc=x(tri(:,3),:);
    

    [dG1,dG2,dG3,dD1,dD2,dD3] = dcross((Xb-Xa)/2,(Xc-Xa)/2,DXig);
    
    dxg2 = U2 + [accumarray(tri(:),[-dG1-dD1;+dG1 ;+dD1 ],[nx,1],[],0),...
                 accumarray(tri(:),[-dG2-dD2;+dG2 ;+dD2 ],[nx,1],[],0),...
                 accumarray(tri(:),[-dG3-dD3;+dG3 ;+dD3 ],[nx,1],[],0)];
    
end

end

