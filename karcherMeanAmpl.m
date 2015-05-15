function [C,vel,grad_norm]=karcherMeanAmpl(d,splineData,quadData,quadDataTensor)
% Compute the Karcher mean of the curves d1,...dn
% Call function as
% Karcher_ampl(d,splineData,quadData,quadDataTensor)
% Input:
%       d, cell array of curves
%       ( d{i}, [NxdSpace], first set of control points)
%       splineData,
%       quadData,
%       quadDataTensor
%
%       
% Output:
%       E_geo, optimal energy
%       dPath, optimal path
%

%Read out number of curves splineData and quadData
Nsteps=50;
n = length(d);
sz = size(d{1});
noCurves = length(d);
vel_mean = zeros(sz);
stop=0;

%Initial guess for Karcher mean
C = d{1}; 

%% Write datfile2 that is valid for the whole minimization
datfile2='H2_tensor.dat';
disp(['main.m, calling writeDatFile2.m, datfile = ' datfile2]);
tic
writeDatFile2(splineData,quadData,quadDataTensor,datfile2);
toc
j=3;

grad_stop = 1e-1;

iter = 0;
max_iter = 20;
while (stop == 0 && iter < max_iter);
    iter = iter + 1;
    %% Calculate the geodesics from C to all curves
    for i=1:n;
        disp(['Iter: ',num2str(iter),' Curve: ',num2str(i)]);
        [~, dPath] = geodesicBvpAmpl(C,d{i},splineData,quadData,quadDataTensor,'datfileexists',true);
        vel_i = pathVelocity(dPath,0,splineData);
        size(vel_i)
        size(vel_mean)
        vel_mean = vel_mean + 1/n*vel_i;
    end
    norm_vel = norm(vel_mean,2);
    disp(['Iter: ',num2str(iter),' Calling geodesicForward']);
    disp(['||grad f|| is ',num2str(norm_vel)]);
    
    if norm_vel < grad_stop;
        stop = 1;
        vel = vel_mean;
        grad_norm = norm_vel;
    else
        C = geodesicForward(C,C + 1/Nsteps*vel_mean,Nsteps+1,...
            splineData,quadData,'endpoint');
        vel_mean = zeros(sz);
    end   
    
end
end

