function [C,vel,norm]=Karcher_ampl(varargin)
% Compute the Karcher mean of the curves d1,...dn
% Call function as
% Karcher_ampl(d1,...,dn,splineData,quadData,quadDataTensor)
% Input:
%       d0, [NxdSpace], first set of control points
%       .
%       .
%       .
%       dn, [NxdSpace], n-th set of control points
%
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
Nsteps=10;
n = nargin-3;
splineData = varargin{n+1};
quadData = varargin{n+2};
quadDataTensor = varargin{n+3};
sz = size(varargin{1});
vel_mean = zeros(sz);
stop=0;

%Initial guess for Karcher mean
C = varargin{1};

%% Write datfile2 that is valid for the whole minimization
datfile2='H2_tensor.dat';
disp(['main.m, calling writedatfile2.m, datfile = ' datfile2]);
tic
writedatfile2(splineData,quadData,quadDataTensor,datfile2);
toc
j=3;


while stop == 0;
%% Calculate the geodesics from C to all curves
    for i=1:n;
    [~, dPath] = geodesicBVP_ampl(C,varargin{i},splineData,quadData,quadDataTensor,0);
    vel_i = pathVelocity(dPath,0, splineData )
    vel_mean = vel_mean + 1/n*vel_i;
    end
    norm_vel = j-1; 
    %norm(vel_mean);
    if norm_vel < 1;
        stop = 1;
        vel = vel_mean;
        norm = norm_vel;
    else
        C = geodesicForward(C,vel_mean,Nsteps);
        vel_mean = zeros(sz);
    end    
end
end

