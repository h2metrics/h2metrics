function [dxg]= dvarifoldnorm(fs1,fs2,objfun)



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
               
        scp_unit_normalsXX = zeros(Tx);
        scp_unit_normalsXY = zeros(Tx,Ty);
    
        for l=1:d
            distance2_center_faceXX = distance2_center_faceXX+(repmat(center_faceX(:,l),1,Tx)-repmat(center_faceX(:,l)',Tx,1)).^2;
            distance2_center_faceXY = distance2_center_faceXY+(repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)).^2;

            scp_unit_normalsXX = scp_unit_normalsXX + (repmat(unit_normalsX(:,l),1,Tx).*repmat(unit_normalsX(:,l)',Tx,1));
            scp_unit_normalsXY = scp_unit_normalsXY + (repmat(unit_normalsX(:,l),1,Ty).*repmat(unit_normalsY(:,l)',Tx,1));
        end

% compute  Geometric kernel      
        Kernel_geomXX  = radial_function_geom(distance2_center_faceXX,0,objfun);
        Kernel_geomXY  = radial_function_geom(distance2_center_faceXY,0,objfun);
        dKernel_geomXX = radial_function_geom(distance2_center_faceXX,1,objfun);
        dKernel_geomXY = radial_function_geom(distance2_center_faceXY,1,objfun);                
        
% compute tangent space kernel
        Kernel_tanXX  = radial_function_sphere(scp_unit_normalsXX,0,objfun);
        Kernel_tanXY  = radial_function_sphere(scp_unit_normalsXY,0,objfun);
        dKernel_tanXX = radial_function_sphere(scp_unit_normalsXX,1,objfun);
        dKernel_tanXY = radial_function_sphere(scp_unit_normalsXY,1,objfun);

% compute Area 
        AreaXX = (norm_normalsX * norm_normalsX');
        AreaXY = (norm_normalsX * norm_normalsY');

%------------------------------------
% Compute derivative wrt center_face 
%------------------------------------
        DXX  = AreaXX .* dKernel_geomXX .* Kernel_tanXX;
        DXY  = AreaXY .* dKernel_geomXY .* Kernel_tanXY;
        
        Dcenter_faceXg=zeros(Tx,d);   
        for l=1:d
            Dcenter_faceXg(:,l)=4*(sum(( DXX .*(repmat(center_faceX(:,l),1,Tx)-repmat(center_faceX(:,l)',Tx,1))),2)...
                -sum(DXY .* (repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)),2)); % scalar kernel case
        end
        
%--------------------------------
% Compute derivative wrt normals 
%--------------------------------
        MXX = Kernel_geomXX .* dKernel_tanXX;
        MXY = Kernel_geomXY .* dKernel_tanXY;                    
        mXX = Kernel_geomXX .* Kernel_tanXX;
        mXY = Kernel_geomXY .* Kernel_tanXY;                
         
        DnormalsXg = 2*(  repmat(mXX * norm_normalsX  - (MXX .* scp_unit_normalsXX  ) * norm_normalsX,1,d) .* unit_normalsX ...
                	+ MXX * normalsX ...
                	- repmat(mXY * norm_normalsY  - (MXY .* scp_unit_normalsXY  ) * norm_normalsY,1,d) .* unit_normalsX...
                	- MXY * normalsY);
            
%end
        
%------------------------        
% End of the chain's rule
%------------------------
        [dxg]= chain_rule(fs1.x,fs1.G,Dcenter_faceXg,DnormalsXg);
         
end
