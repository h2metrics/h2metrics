function [g, dxg1, dxg2]= dvarifoldnorm(fs1, fs2, objfun)

d=size(fs1.x,2);

[center_faceX,normalsX]=varatoms(fs1.x,fs1.G);
[center_faceY,normalsY]=varatoms(fs2.x,fs2.G);

% Compute unit normals
Tx=size(center_faceX,1);
Ty=size(center_faceY,1);

% Compute unit normals
norm_normalsX = sqrt(sum(normalsX .^2,2));
norm_normalsY = sqrt(sum(normalsY .^2,2));

unit_normalsX = normalsX ./  repmat(norm_normalsX,1,size(normalsX,2));
unit_normalsY = normalsY ./  repmat(norm_normalsY,1,size(normalsY,2));

% switch objfun.method
%     case 'cuda'  % use cuda to speedup the computation
% 
% 	 eval(['dffshape_scp=@dffshape_dist_Gpu_',lower(objfun.kernel_geom),lower(objfun.kernel_signal),lower(objfun.kernel_grass),';']);
%         
%         DsignalXg = 2 *(dffshape_scp(center_faceX',center_faceX',signalX',signalX',normalsX',normalsX',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass)...
%             - dffshape_scp(center_faceX',center_faceY',signalX',signalY',normalsX',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass))';
%         
% 
% 	 eval(['dXifshape_scp=@dXifshape_dist_Gpu_',lower(objfun.kernel_geom),lower(objfun.kernel_signal),lower(objfun.kernel_grass),';']);
%         DnormalsXg = 2 * (dXifshape_scp(center_faceX',center_faceX',signalX',signalX',normalsX',normalsX',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass)...
%             - dXifshape_scp(center_faceX',center_faceY',signalX',signalY',normalsX',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass))';
%         
% 	 eval(['dXfshape_scp=@dXfshape_dist_Gpu_',lower(objfun.kernel_geom),lower(objfun.kernel_signal),lower(objfun.kernel_grass),';']);
%         Dcenter_faceXg = 2 * (dXfshape_scp(center_faceX',center_faceX',signalX',signalX',normalsX',normalsX',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass)...
%             - dXfshape_scp(center_faceX',center_faceY',signalX',signalY',normalsX',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass))';
%         
%     otherwise


% compute squared distances and angles
 
distance2_center_faceXX = zeros(Tx);
distance2_center_faceXY=zeros(Tx,Ty);
distance2_center_faceYY=zeros(Ty);

scp_unit_normalsXX = zeros(Tx);
scp_unit_normalsXY = zeros(Tx,Ty);
scp_unit_normalsYY = zeros(Ty);

for l=1:d        
    distance2_center_faceXX = distance2_center_faceXX+(repmat(center_faceX(:,l),1,Tx)-repmat(center_faceX(:,l)',Tx,1)).^2;
    distance2_center_faceXY = distance2_center_faceXY+(repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)).^2;
    distance2_center_faceYY = distance2_center_faceYY+(repmat(center_faceY(:,l),1,Ty)-repmat(center_faceY(:,l)',Ty,1)).^2;

    scp_unit_normalsXX = scp_unit_normalsXX + (repmat(unit_normalsX(:,l),1,Tx).*repmat(unit_normalsX(:,l)',Tx,1));
    scp_unit_normalsXY = scp_unit_normalsXY + (repmat(unit_normalsX(:,l),1,Ty).*repmat(unit_normalsY(:,l)',Tx,1));
    scp_unit_normalsYY = scp_unit_normalsYY + (repmat(unit_normalsY(:,l),1,Ty).*repmat(unit_normalsY(:,l)',Ty,1));
end
scp_unit_normalsYX = scp_unit_normalsXY';
        
% compute  Geometric kernel      
[Kernel_geomXX, dKernel_geomXX] = ...
    radial_function_geom(distance2_center_faceXX,objfun);
[Kernel_geomXY, dKernel_geomXY] = ...
    radial_function_geom(distance2_center_faceXY,objfun);

Kernel_geomYX  = Kernel_geomXY';
dKernel_geomYX = dKernel_geomXY';

[Kernel_geomYY, dKernel_geomYY] = ...
    radial_function_geom(distance2_center_faceYY,objfun);

% compute tangent space kernel
[Kernel_tanXX, dKernel_tanXX] = ...
    radial_function_sphere(scp_unit_normalsXX,objfun);
[Kernel_tanXY, dKernel_tanXY] = ...
    radial_function_sphere(scp_unit_normalsXY,objfun);
Kernel_tanYX = Kernel_tanXY';
dKernel_tanYX = dKernel_tanXY';
[Kernel_tanYY, dKernel_tanYY] = ...
    radial_function_sphere(scp_unit_normalsYY,objfun);

% compute Area 
AreaXX = (norm_normalsX * norm_normalsX');
AreaXY = (norm_normalsX * norm_normalsY');
AreaYX = AreaXY';
AreaYY = (norm_normalsY * norm_normalsY');

% norm(x)=
PXX = sum(sum(AreaXX .* Kernel_geomXX .* Kernel_tanXX ));

% morm(y) =
PYY =sum(sum(AreaYY .* Kernel_geomYY .* Kernel_tanYY));

%prs(x,y) =
PXY = sum(sum(AreaXY .* Kernel_geomXY .* Kernel_tanXY));

% Varifold distance
g = PXX + PYY - 2* PXY;

%------------------------------------
% Compute derivative wrt center_face 
%------------------------------------
DXX  = AreaXX .* dKernel_geomXX .* Kernel_tanXX;
DXY  = AreaXY .* dKernel_geomXY .* Kernel_tanXY;
DYX  = AreaYX .* dKernel_geomYX .* Kernel_tanYX;
DYY  = AreaYY .* dKernel_geomYY .* Kernel_tanYY;

Dcenter_faceXg=zeros(Tx,d);   
Dcenter_faceYg=zeros(Ty,d);
for l=1:d            
    Dcenter_faceXg(:,l)=4*(sum(( DXX .*(repmat(center_faceX(:,l),1,Tx)-repmat(center_faceX(:,l)',Tx,1))),2)...
        -sum(DXY .* (repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)),2)); % scalar kernel case
    Dcenter_faceYg(:,l)=4*(sum(( DYY .*(repmat(center_faceY(:,l),1,Ty)-repmat(center_faceY(:,l)',Ty,1))),2)...
        -sum(DYX .* (repmat(center_faceY(:,l),1,Tx)-repmat(center_faceX(:,l)',Ty,1)),2)); % scalar kernel case
end
        
%--------------------------------
% Compute derivative wrt normals 
%--------------------------------
MXX = Kernel_geomXX .* dKernel_tanXX;
MXY = Kernel_geomXY .* dKernel_tanXY;                    
MYX = Kernel_geomYX .* dKernel_tanYX;
MYY = Kernel_geomYY .* dKernel_tanYY;
mXX = Kernel_geomXX .* Kernel_tanXX;
mXY = Kernel_geomXY .* Kernel_tanXY;
mYX = Kernel_geomYX .* Kernel_tanYX;
mYY = Kernel_geomYY .* Kernel_tanYY;

tmpXX = repmat(mXX * norm_normalsX - (MXX .* scp_unit_normalsXX) * norm_normalsX, 1, d) .* unit_normalsX ...
            + MXX * normalsX;
tmpYY = repmat(mYY * norm_normalsY - (MYY .* scp_unit_normalsYY) * norm_normalsY, 1, d) .* unit_normalsY ...
            + MYY * normalsY;
tmpXY = repmat(mXY * norm_normalsY - (MXY .* scp_unit_normalsXY) * norm_normalsY,1,d) .* unit_normalsX...
            + MXY * normalsY;
tmpYX = repmat(mYX * norm_normalsX - (MYX .* scp_unit_normalsYX) * norm_normalsX,1,d) .* unit_normalsY...
            + MYX * normalsX;

DnormalsXg = 2 * (tmpXX - tmpXY);
DnormalsYg = 2 * (tmpYY - tmpYX);
        
% DnormalsXg = 2*(  repmat(mXX * norm_normalsX  - (MXX .* scp_unit_normalsXX  ) * norm_normalsX,1,d) .* unit_normalsX ...
%             + MXX * normalsX ...
%             - repmat(mXY * norm_normalsY  - (MXY .* scp_unit_normalsXY  ) * norm_normalsY,1,d) .* unit_normalsX...
%             - MXY * normalsY);
% 
% DnormalsYg = 2*(  repmat(mYY * norm_normalsY  - (MYY .* scp_unit_normalsYY  ) * norm_normalsY,1,d) .* unit_normalsY ...
%             + MYY * normalsY ...
%             - repmat(mYX * norm_normalsX  - (MYX .* scp_unit_normalsYX  ) * norm_normalsX,1,d) .* unit_normalsY...
%             - MYX * normalsX);

%----------------------
% End of the chain rule
%----------------------
[dxg1]= chain_rule(fs1.x,fs1.G,Dcenter_faceXg,DnormalsXg);
[dxg2]= chain_rule(fs2.x,fs2.G,Dcenter_faceYg,DnormalsYg);
         
end
