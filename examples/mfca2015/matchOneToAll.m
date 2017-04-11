function [Dist,InitialVel] = matchOneToAll(d0,d,splineData,quadData,quadDataTensor)
% Matches one curve to all curves d1,...dn
% Call function as
% Karcher_ampl(d,splineData,quadData,quadDataTensor)
% Input:
%       d0, curve
%       d, cell array of curves
%       ( d{i}, [NxdSpace], first set of control points)
%       splineData,
%       quadData,
%       quadDataTensor
%
%       
% Output:
%       distanceVector

%Read out number of curves splineData and quadData
n = length(d);
Dist = zeros(n,1);

%Initial guess for Karcher mean

%% Write datfile2 that is valid for the whole minimization
datfile2='H2_tensor.dat';
disp(['main.m, calling writeDatFile2.m, datfile = ' datfile2]);
tic
writeDatFile2(splineData,quadData,quadDataTensor,datfile2);
toc

for i=1:n;
        disp([' Curve: ',num2str(i)]);
        [~,dPath] = geodesicBvpAmpl(d0,d{i},splineData,quadData,quadDataTensor,'datfileexists',true);
        InitialVel{i} = pathVelocity(dPath,0,splineData);
        Dist(i,1) = pathRiemH2Length(dPath,splineData,quadData,quadDataTensor);
end  
end

