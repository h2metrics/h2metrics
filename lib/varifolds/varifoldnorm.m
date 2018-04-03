function g=varifoldnorm(fs1,fs2,objfun)
% FVARIFOLDNORM(templatefinal,target,objfun) computes kernel based
% distances between fshapes.
%
%  \sum_i\sum_j K_geom(-norm(x_i-y_j)^2) K_tan(angle(V_i,W_j))
% 
% Possible method are 'cuda' or 'matlab'.
%
% Inputs:
%   fs1 : structure containing the first shape (source)
%   fs2 : structure containing the target.
%   objfun : is a structure containing the data attachment term parameters
% Output
%   g : a real number.
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)

d=size(fs1.x,2);

% discretize the fshapes
[center_faceX,normalsX]=varatoms(fs1.x,fs1.G);
[center_faceY,normalsY]=varatoms(fs2.x,fs2.G);

%switch objfun.method
%      case 'cuda'  % use cuda to speedup the computation
% 
% 	 eval(['fshape_scp=@fshape_dist_Gpu_',lower(objfun.kernel_geom),lower(objfun.kernel_signal),lower(objfun.kernel_grass),';']);
% 
%          % norm(x)^2 =
%          XX= fshape_scp(center_faceX',center_faceX',signalX',signalX',normalsX',normalsX',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass);
%          PXX =  sum(XX);
%          
%          
%          % morm(y)^2 =
%          YY= fshape_scp(center_faceY',center_faceY',signalY',signalY',normalsY',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass);
%          PYY =  sum(YY);
%          
%          %prs(x,y) =
%          XY= fshape_scp(center_faceX',center_faceY',signalX',signalY',normalsX',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass);
%          PXY =  sum(XY);
%                   
%      otherwise

% Compute unit normals

        Nx=size(center_faceX,1);
        Ny=size(center_faceY,1);
        % Compute unit normals
        norm_normalsX = sqrt(sum(normalsX .^2,2));
        norm_normalsY = sqrt(sum(normalsY .^2,2));
        
        unit_normalsX = normalsX ./  repmat(norm_normalsX,1,size(normalsX,2));
        unit_normalsY = normalsY ./  repmat(norm_normalsY,1,size(normalsY,2));
        
        %compute squared distances and angles       
        distance_center_faceXX = zeros(Nx);
        distance_center_faceYY=zeros(Ny);
        distance_center_faceXY=zeros(Nx,Ny);
                
        scp_unit_normalsXX = zeros(Nx);
        scp_unit_normalsYY = zeros(Ny);        
        scp_unit_normalsXY = zeros(Nx,Ny);
    
        for l=1:d
            distance_center_faceXX = distance_center_faceXX+(repmat(center_faceX(:,l),1,Nx)-repmat(center_faceX(:,l)',Nx,1)).^2;
            distance_center_faceYY = distance_center_faceYY+(repmat(center_faceY(:,l),1,Ny)-repmat(center_faceY(:,l)',Ny,1)).^2;
            distance_center_faceXY = distance_center_faceXY+(repmat(center_faceX(:,l),1,Ny)-repmat(center_faceY(:,l)',Nx,1)).^2;

            scp_unit_normalsXX = scp_unit_normalsXX + (repmat(unit_normalsX(:,l),1,Nx).*repmat(unit_normalsX(:,l)',Nx,1));
            scp_unit_normalsYY = scp_unit_normalsYY + (repmat(unit_normalsY(:,l),1,Ny).*repmat(unit_normalsY(:,l)',Ny,1));
            scp_unit_normalsXY = scp_unit_normalsXY + (repmat(unit_normalsX(:,l),1,Ny).*repmat(unit_normalsY(:,l)',Nx,1));
        end
        
        % Geometric kernel      
        Kernel_geomXX = radial_function_geom(distance_center_faceXX,objfun);
        Kernel_geomYY = radial_function_geom(distance_center_faceYY,objfun);
        Kernel_geomXY = radial_function_geom(distance_center_faceXY,objfun);
        

        % tangent space kernel
        Kernel_tanXX = radial_function_sphere(scp_unit_normalsXX,objfun);
        Kernel_tanYY =  radial_function_sphere(scp_unit_normalsYY,objfun);
        Kernel_tanXY = radial_function_sphere(scp_unit_normalsXY,objfun);
        
        % Area 
        AreaXX = (norm_normalsX * norm_normalsX');
        AreaYY = (norm_normalsY * norm_normalsY');
        AreaXY = (norm_normalsX * norm_normalsY');
        
        % norm(x)=
        PXX = sum(sum(AreaXX .* Kernel_geomXX .* Kernel_tanXX ));
        
        % morm(y) =
        PYY =sum(sum(AreaYY .* Kernel_geomYY .* Kernel_tanYY));
        
        %prs(x,y) =
        PXY = sum(sum(AreaXY .* Kernel_geomXY .* Kernel_tanXY));
        
%end         

 g= PXX + PYY - 2* PXY;
end

%------------------------------------%
% equivalent code with scalarProduct %
%------------------------------------%


%function g = shape_Kernel_distance(fs1,fs2,objfun)
%% 
%% Possible method are 'cuda' or 'matlab'.

%% discretize the fshapes
%[center_faceX,signalX,normalsX]=fcatoms(fs1.x,fs1.f,fs1.G,objfun.signal_type);
%[center_faceY,signalY,normalsY]=fcatoms(fs2.x,fs2.f,fs2.G,objfun.data_signal_type);

%% compute the morm squared
%PXX = shape_Kernel_scalarProduct(center_faceX,signalX,normalsX,center_faceX,signalX,normalsX,objfun);
%PYY = shape_Kernel_scalarProduct(center_faceY,signalY,normalsY,center_faceY,signalY,normalsY,objfun);
%PXY = shape_Kernel_scalarProduct(center_faceX,signalX,normalsX,center_faceY,signalY,normalsY,objfun);

%g= PXX + PYY - 2* PXY;

%end




%function res = shape_Kernel_scalarProduct(center_faceX,signalX,normalsX,center_faceY,signalY,normalsY,objfun)
%% SHAPE_KERNEL_SCALARPRODUCT(templatefinal,target,objfun) computes the functional 
%% varifold distance with Gaussian kernel :
%%
%%  \sum_i\sum_j K_signal(f_i-g_j)^2) K_geom(-norm(x_i-y_j)^2) K_tan(angle(V_i,W_j))
%%
%% Inputs:
%%   center_faceX : matrix Tx * d, containing the coordinates of the center of faces.
%%   normalsX : matrix Tx * d, containing the coordinates of the normal to the faces.
%%   signalX : matrix Tx *1 (column vector) containing the signals.
%%   center_faceY : matrix Ty * d, containing the coordinates of the center of faces.
%%   normalsY : matrix Ty * d, containing the coordinates of the normal to the faces.
%%   signalY : matrix Ty *1 (column vector) containing the signals.
%%   objfun : structure containing the options for the kernels (geom, signal and tan)
%% Output
%%   g : a real number.
%%
%% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, G. Nardi, A. Trouve (2012-2014)

%% Get dimensions
%d=size(center_faceX,2);
%Tx=size(center_faceX,1);
%Ty=size(center_faceY,1);

%% Compute norms of the normals
%norm_normalsX = sqrt(sum(normalsX .^2,2));
%norm_normalsY = sqrt(sum(normalsY .^2,2));

%% Compute unit normals
%unit_normalsX = normalsX ./  repmat(norm_normalsX,1,size(normalsX,2));
%unit_normalsY = normalsY ./  repmat(norm_normalsY,1,size(normalsY,2));

%%compute squared distances and angles
%distance_signalXY = (repmat(signalX,1,Ty)-repmat(signalY',Tx,1)).^2;        
%distance_center_faceXY=zeros(Tx,Ty);
%oriented_angle_normalsXY = zeros(Tx,Ty);

%for l=1:d
    %distance_center_faceXY = distance_center_faceXY+(repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)).^2;
    %oriented_angle_normalsXY = oriented_angle_normalsXY + (repmat(unit_normalsX(:,l),1,Ty).*repmat(unit_normalsY(:,l)',Tx,1));
%end

%% Kernels
%Kernel_geomXY = radial_function_geom(distance_center_faceXY,0,objfun);
%Kernel_signalXY = radial_function_signal(distance_signalXY,0,objfun);
%Kernel_tanXY = radial_function_sphere(oriented_angle_normalsXY,0,objfun);

%%prs(x,y) =
%res = sum(sum((norm_normalsX * norm_normalsY') .* Kernel_geomXY .* Kernel_signalXY .* Kernel_tanXY));

%end



