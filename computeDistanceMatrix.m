function Dist=computeDistanceMatrix(d,splineData,quadData,quadDataTensor,varargin)
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
%       MAtrix of pairwise distances
%

%Read out number of curves splineData and quadData
n = length(d);
Dist = zeros(n);


%% Write datfile2 that is valid for the whole minimization
datfile2='H2_tensor.dat';
disp(['main.m, calling writeDatFile2.m, datfile = ' datfile2]);
tic
writeDatFile2(splineData,quadData,quadDataTensor,datfile2);
toc


for i=1:(n-1);
    for j=(i+1):n;
     [~, dPath] = geodesicBvpAmpl(d{i},d{j},splineData,quadData,quadDataTensor,'datfileexists',true);
     Dist(i,j) = energyH2(dPath,splineData,quadDataTensor);
     Dist(j,i) = Dist(i,j);
    end
end 
end

