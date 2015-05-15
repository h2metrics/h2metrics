function [L2coeff H1coeff H2coeff]=determineConstants(d,splineData,quadData,quadDataTensor)
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
L2value=0;
H1value=0;
H2value=0;

for i=1:(n-1);
    for j=(i+1):n;
    dPath = linearPath(d{i},d{j},splineData);
    [~,Value]=pathRiemH2Length( dPath,splineData,quadData, quadDataTensor);
    L2value=L2value+Value(1);
    H1value=H1value+Value(2);
    H2value=H2value+Value(3);
    end
end 

L2coeff = 1;
H1coeff = L2value/H1value;
H2coeff = L2value/H2value;
end

