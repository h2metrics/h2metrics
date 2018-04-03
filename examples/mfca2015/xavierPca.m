
%% PCA using all Patients from allpatients mean
clear all;
load('xavier_splines_v1_N30_noloops','dPatients_noloops');
splineData = constructEmptySplineData;
splineData.N = 30; %no. control points, must be bigger than n+1
splineData.Nt = 10 + 2; %Number of time control points
splineData.Nphi = 6; %No. control points for diffeomorphisms
splineData.nS = 3; %spacial degree
splineData.nT = 2; %time degree
splineData.nPhi = 3; %diffemorphism degree
splineData.quadDegree = [8,4]; %Quadrature precission
splineData.dSpace = 2;
splineData.noInterpolS = 2 * splineData.N; % For composition
splineData = constructKnots(splineData);
[quadData, quadDataTensor] = setupQuadData(splineData);
dSick = {dPatients_noloops{10:19}};
dHealthy = {dPatients_noloops{1:9}};
splineData.a = [1 1 1];
noCurves=19;
noAllCurves=21;
dSpace=2;
n1 = 9;
n2 = 10;
N=splineData.N;
load('InitialVelPatientsA1B1C1.mat');
load('CPatientsA1B1C1.mat');
V = zeros(splineData.N*splineData.dSpace,noCurves);
for ii = 1:noCurves;
    V(:,ii) = InitialVelPatientsA1B1C1{ii}(:);
end
G = metricMatrixH2(CPatients,splineData,quadData);
Gsqrt = sqrtm(G);
Sigma = zeros([N*dSpace, N*dSpace]);
for jj = 1:noCurves;
    Sigma = Sigma + V(:,jj) * V(:,jj)';
end


Sigma = 1./(noCurves-1) * Gsqrt * Sigma * Gsqrt;

[U, Lambda] = eig(Sigma);
Lambda = real(diag(Lambda))

% Calculate 2d coordinates (also of the means)
for ii = 1:19;
    V(:,ii) = InitialVelPatientsA1B1C1{ii}(:);
end
W = U(:,1:2);
pts2d = zeros([2, 21]);
Vproj = zeros([2, 21]);
sz = size(V(:,1));

for jj = 2:n1+1
    pts2d(:,jj) = W' * Gsqrt * V(:,jj-1);
    Vproj(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
end
for jj = n1+2:n1+n2+1
    pts2d(:,jj) = W' * Gsqrt * V(:,jj-1);
    Vproj(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
end


%pts2d(:,1) = W' * Gsqrt * VHealthy;
%Vproj(:,1) = pts2d(:,1) ./ sqrt(Lambda(1:2));
%pts2d(:,noAllCurves) = W' * Gsqrt * VSick;
%Vproj(:,noAllCurves) = pts2d(:,noAllCurves) ./ sqrt(Lambda(1:2));


%Plot
figure(2)
hold on;
axis square;
axis([-1.5 1.5 -2.5 3.5])
grid on;
box on;
Y = Vproj';
%plot(0,0,'o')
%text(0, 0, ['mP'],'Color','black','verticalAlignment', 'middle');

plot(Y(2:n1+1,1),Y(2:n1+1,2), 'o','Color','blue');
plot(Y(n1+2:n1+n2+1,1),Y(n1+2:n1+n2+1,2), 'x','Color','red');

%plot(Y(1,1),Y(1,2), 'o','Color','blue');
plot(0,0, 'x','Color','black');
for kk = 2:n1+1
    text(Y(kk,1), Y(kk,2), [' ',num2str(kk-1)],'Color','blue','verticalAlignment', 'middle');
end
for kk = n1+2:n1+n2+1
    text(Y(kk,1), Y(kk,2), [' ',num2str(kk-1)],'Color','red','verticalAlignment', 'middle');
end

%text(Y(1,1), Y(1,2), ' H','Color','blue','verticalAlignment', 'middle');
text(0, 0, ' M','Color','black','verticalAlignment', 'middle');
pbaspect();
daspect();
set(gcf,'color','w');
hold off;
export_fig xavierPca.pdf




%% PCA only for the Healthy Patients
clear all;
load('xavier_splines_v1_N30_noloops','dPatients_noloops');
splineData = constructEmptySplineData;
splineData.N = 30; %no. control points, must be bigger than n+1
splineData.Nt = 10 + 2; %Number of time control points
splineData.Nphi = 6; %No. control points for diffeomorphisms
splineData.nS = 3; %spacial degree
splineData.nT = 2; %time degree
splineData.nPhi = 3; %diffemorphism degree
splineData.quadDegree = [8,4]; %Quadrature precission
splineData.dSpace = 2;
splineData.noInterpolS = 2 * splineData.N; % For composition
splineData = constructKnots(splineData);
[quadData, quadDataTensor] = setupQuadData(splineData);
dSick = {dPatients_noloops{10:19}};
dHealthy = {dPatients_noloops{1:9}};
splineData.a = [1 1 1];

close all;
noCurves=9;
n1=9;
dSpace=2;
N=splineData.N;
load('InitialVelHealthy.mat');
load('CHealthyA1B1C1.mat');
V = zeros(splineData.N*splineData.dSpace,noCurves);
for ii = 1:noCurves;
    V(:,ii) = InitialVelHealthy{ii}(:);
end
G = metricMatrixH2(CHealthy,splineData,quadData);
Gsqrt = sqrtm(G);
Sigma = zeros([N*dSpace, N*dSpace]);
for jj = 1:noCurves;
    Sigma = Sigma + V(:,jj) * V(:,jj)';
end


Sigma = 1./(noCurves-1) * Gsqrt * Sigma * Gsqrt;

[U, Lambda] = eig(Sigma);
Lambda = real(diag(Lambda))

% Calculate 2d coordinates (also of the means)
for ii = 1:9;
    V(:,ii) = InitialVelHealthy{ii}(:);
end
W = U(:,1:2);
pts2d = zeros([2, 21]);
Vproj = zeros([2, 21]);
sz = size(V(:,1));
VHealthy=zeros(sz);
for jj = 2:10
    pts2d(:,jj) = W' * Gsqrt * V(:,jj-1);
    Vproj(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
end


pts2d(:,1) = W' * Gsqrt * VHealthy;
Vproj(:,1) = pts2d(:,1) ./ sqrt(Lambda(1:2));


%Plot
hold on;
axis square;
axis([-2 2 -2.5 1.5])
grid on;
box on;
Y = Vproj';
%plot(0,0,'o')
%text(0, 0, ['mP'],'Color','black','verticalAlignment', 'middle');

plot(Y(2:n1+1,1),Y(2:n1+1,2), 'o','Color','blue');

plot(Y(1,1),Y(1,2), 'o','Color','black');

for kk = 2:n1+1
    text(Y(kk,1), Y(kk,2), [' ',num2str(kk-1)],'Color','blue','verticalAlignment', 'middle');
end

text(Y(1,1), Y(1,2), ' H','Color','black','verticalAlignment', 'middle');
set(gcf,'color','w');
hold off;
export_fig xavierPcaHealthyOnly.pdf


%% PCA only for Sick Patients
clear all;
load('xavier_splines_v1_N30_noloops','dPatients_noloops');
splineData = constructEmptySplineData;
splineData.N = 30; %no. control points, must be bigger than n+1
splineData.Nt = 10 + 2; %Number of time control points
splineData.Nphi = 6; %No. control points for diffeomorphisms
splineData.nS = 3; %spacial degree
splineData.nT = 2; %time degree
splineData.nPhi = 3; %diffemorphism degree
splineData.quadDegree = [8,4]; %Quadrature precission
splineData.dSpace = 2;
splineData.noInterpolS = 2 * splineData.N; % For composition
splineData = constructKnots(splineData);
[quadData, quadDataTensor] = setupQuadData(splineData);
dSick = {dPatients_noloops{10:19}};
dHealthy = {dPatients_noloops{1:9}};
splineData.a = [1 1 1];

close all;
noCurves=10;
dSpace=2;
n1 = 9;
n2 = 10;
N=splineData.N;
load('InitialVelSick.mat');
load('CSickA1B1C1.mat');
V = zeros(splineData.N*splineData.dSpace,noCurves);
for ii = 1:noCurves;
    V(:,ii) = InitialVelSick{9+ii}(:);
end
G = metricMatrixH2(CSick,splineData,quadData);
Gsqrt = sqrtm(G);
Sigma = zeros([N*dSpace, N*dSpace]);
for jj = 1:noCurves;
    Sigma = Sigma + V(:,jj) * V(:,jj)';
end


Sigma = 1./(noCurves-1) * Gsqrt * Sigma * Gsqrt;

[U, Lambda] = eig(Sigma);
Lambda = real(diag(Lambda));

% Calculate 2d coordinates (also of the means)
for ii = 1:19;
    V(:,ii) = InitialVelSick{ii}(:);
end
W = U(:,1:2);
pts2d = zeros([2, 21]);
Vproj = zeros([2, 21]);
sz = size(V(:,1));
VSick=zeros(sz);

for jj = n1+2:n1+n2+1
    pts2d(:,jj) = W' * Gsqrt * V(:,jj-1);
    Vproj(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
end


pts2d(:,21) = W' * Gsqrt * VSick;
Vproj(:,21) = pts2d(:,noCurves+1) ./ sqrt(Lambda(1:2));


%Plot
hold on;
axis square;
axis([-2.5 1.5 -2.25 1.75])
grid on;
box on;
Y = Vproj';
%plot(0,0,'o')
%text(0, 0, ['mP'],'Color','black','verticalAlignment', 'middle');

plot(Y(n1+2:n1+n2+1,1),Y(n1+2:n1+n2+1,2), 'x','Color','red');
plot(0,0, 'x','Color','black');

for kk = n1+2:n1+n2+1
    text(Y(kk,1), Y(kk,2), [' ',num2str(kk-1)],'Color','red','verticalAlignment', 'middle');
end
text(0, 0, ' S','Color','black','verticalAlignment', 'middle');

set(gcf,'color','w');
hold off;
export_fig xavierPcaSickOnly.pdf






% %% PCA using Healthy Patients
% clear all;
% load('xavier_splines_v1_N30_noloops','dPatients_noloops');
% splineData = constructEmptySplineData;
% splineData.N = 30; %no. control points, must be bigger than n+1
% splineData.Nt = 10 + 2; %Number of time control points
% splineData.Nphi = 6; %No. control points for diffeomorphisms
% splineData.nS = 3; %spacial degree
% splineData.nT = 2; %time degree
% splineData.nPhi = 3; %diffemorphism degree
% splineData.quadDegree = [8,4]; %Quadrature precission
% splineData.dSpace = 2;
% splineData.noInterpolS = 2 * splineData.N; % For composition
% splineData = constructKnots(splineData);
% [quadData, quadDataTensor] = setupQuadData(splineData);
% dSick = {dPatients_noloops{10:19}};
% dHealthy = {dPatients_noloops{1:9}};
% splineData.a = [1 1 1];
% 
% close all;
% noCurves=9;
% noAllCurves=21;
% dSpace=2;
% n1 = 9;
% n2 = 10;
% N=splineData.N;
% load('InitialVelHealthy.mat');
% load('InitialVelMeanToMean.mat');
% load('CHealthyA1B1C1.mat');
% load('CSickA1B1C1.mat');
% V = zeros(splineData.N*splineData.dSpace,noCurves);
% for ii = 1:noCurves;
%     V(:,ii) = InitialVelHealthy{ii}(:);
% end
% G = metricMatrixH2(CHealthy,splineData,quadData);
% Gsqrt = sqrtm(G);
% Sigma = zeros([N*dSpace, N*dSpace]);
% for jj = 1:noCurves;
%     Sigma = Sigma + V(:,jj) * V(:,jj)';
% end
% 
% 
% Sigma = 1./(noCurves-1) * Gsqrt * Sigma * Gsqrt;
% 
% [U, Lambda] = eig(Sigma);
% Lambda = real(diag(Lambda));
% 
% % Calculate 2d coordinates (also of the means)
% for ii = 1:19;
%     V(:,ii) = InitialVelHealthy{ii}(:);
% end
% W = U(:,1:2);
% pts2d = zeros([2, 21]);
% Vproj = zeros([2, 21]);
% sz = size(V(:,1));
% VHealthy=zeros(sz);
% for jj = 2:n1+1
%     pts2d(:,jj) = W' * Gsqrt * V(:,jj-1);
%     Vproj(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
% end
% for jj = n1+2:n1+n2+1
%     pts2d(:,jj) = W' * Gsqrt * V(:,jj-1);
%     Vproj(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
% end
% VSick = -InitialVelMeanToMean{1}(:);
% 
% pts2d(:,1) = W' * Gsqrt * VHealthy;
% Vproj(:,1) = pts2d(:,1) ./ sqrt(Lambda(1:2));
% pts2d(:,noAllCurves) = W' * Gsqrt * VSick;
% Vproj(:,noAllCurves) = pts2d(:,noAllCurves) ./ sqrt(Lambda(1:2));
% 
% 
% %Plot
% hold on;
% axis square;
% axis([-2 2 -2.5 3.7])
% grid on;
% box on;
% Y = Vproj';
% %plot(0,0,'o')
% %text(0, 0, ['mP'],'Color','black','verticalAlignment', 'middle');
% 
% plot(Y(2:n1+1,1),Y(2:n1+1,2), 'o','Color','blue');
% plot(Y(n1+2:n1+n2+1,1),Y(n1+2:n1+n2+1,2), 'x','Color',[0.5 0.5 0.5]);
% 
% plot(Y(1,1),Y(1,2), 'o','Color','black');
% plot(Y(n1+n2+2,1),Y(n1+n2+2,2), 'x','Color',[0.5 0.5 0.5]);
% for kk = 2:n1+1
%     text(Y(kk,1), Y(kk,2), [' ',num2str(kk-1)],'Color','blue','verticalAlignment', 'middle');
% end
% for kk = n1+2:n1+n2+1
%     text(Y(kk,1), Y(kk,2), [' ',num2str(kk-1)],'Color',[0.5 0.5 0.5],'verticalAlignment', 'middle');
% end
% 
% text(Y(1,1), Y(1,2), ' H','Color','black','verticalAlignment', 'middle');
% text(Y(n1+n2+2,1), Y(n1+n2+2,2), ' S','Color',[0.5 0.5 0.5],'verticalAlignment', 'middle');
% daspect([1.51 1  1])
% set(gcf,'color','w');
% hold off;
% export_fig xavierPcaHealthy.pdf
% 
% %% PCA using Sick Patients
% clear all;
% load('xavier_splines_v1_N30_noloops','dPatients_noloops');
% splineData = constructEmptySplineData;
% splineData.N = 30; %no. control points, must be bigger than n+1
% splineData.Nt = 10 + 2; %Number of time control points
% splineData.Nphi = 6; %No. control points for diffeomorphisms
% splineData.nS = 3; %spacial degree
% splineData.nT = 2; %time degree
% splineData.nPhi = 3; %diffemorphism degree
% splineData.quadDegree = [8,4]; %Quadrature precission
% splineData.dSpace = 2;
% splineData.noInterpolS = 2 * splineData.N; % For composition
% splineData = constructKnots(splineData);
% [quadData, quadDataTensor] = setupQuadData(splineData);
% dSick = {dPatients_noloops{10:19}};
% dHealthy = {dPatients_noloops{1:9}};
% splineData.a = [1 1 1];
% 
% close all;
% noCurves=10;
% noAllCurves=21;
% dSpace=2;
% n1 = 9;
% n2 = 10;
% N=splineData.N;
% load('InitialVelSick.mat');
% load('InitialVelMeanToMean.mat');
% load('CSickA1B1C1.mat');
% V = zeros(splineData.N*splineData.dSpace,noCurves);
% for ii = 1:noCurves;
%     V(:,ii) = InitialVelSick{9+ii}(:);
% end
% G = metricMatrixH2(CSick,splineData,quadData);
% Gsqrt = sqrtm(G);
% Sigma = zeros([N*dSpace, N*dSpace]);
% for jj = 1:noCurves;
%     Sigma = Sigma + V(:,jj) * V(:,jj)';
% end
% 
% 
% Sigma = 1./(noCurves-1) * Gsqrt * Sigma * Gsqrt;
% 
% [U, Lambda] = eig(Sigma);
% Lambda = real(diag(Lambda));
% 
% % Calculate 2d coordinates (also of the means)
% for ii = 1:19;
%     V(:,ii) = InitialVelSick{ii}(:);
% end
% W = U(:,1:2);
% pts2d = zeros([2, 21]);
% Vproj = zeros([2, 21]);
% sz = size(V(:,1));
% VSick=zeros(sz);
% for jj = 2:n1+1
%     pts2d(:,jj) = W' * Gsqrt * V(:,jj-1);
%     Vproj(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
% end
% for jj = n1+2:n1+n2+1
%     pts2d(:,jj) = W' * Gsqrt * V(:,jj-1);
%     Vproj(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
% end
% VHealthy = InitialVelMeanToMean{1}(:);
% 
% pts2d(:,1) = W' * Gsqrt * VHealthy;
% Vproj(:,1) = pts2d(:,1) ./ sqrt(Lambda(1:2));
% pts2d(:,noAllCurves) = W' * Gsqrt * VSick;
% Vproj(:,noAllCurves) = pts2d(:,noAllCurves) ./ sqrt(Lambda(1:2));
% 
% 
% %Plot
% hold on;
% axis equal;
% axis([-2.5 2.5 -10 2])
% aspectRatio = pbaspect;
% grid on;
% box on;
% Y = Vproj';
% %plot(0,0,'o')
% %text(0, 0, ['mP'],'Color','black','verticalAlignment', 'middle');
% 
% plot(Y(2:n1+1,1),Y(2:n1+1,2), 'o','Color',[0.5 0.5 0.5]);
% plot(Y(n1+2:n1+n2+1,1),Y(n1+2:n1+n2+1,2), 'x','Color','red');
% 
% plot(Y(1,1),Y(1,2), 'o','Color',[0.5 0.5 0.5]);
% plot(Y(n1+n2+2,1),Y(n1+n2+2,2), 'x','Color','black');
% for kk = 2:n1+1
%     text(Y(kk,1), Y(kk,2), [' ',num2str(kk-1)],'Color',[0.5 0.5 0.5],'verticalAlignment', 'middle');
% end
% for kk = n1+2:n1+n2+1
%     text(Y(kk,1), Y(kk,2), [' ',num2str(kk-1)],'Color','red','verticalAlignment', 'middle');
% end
% 
% text(Y(1,1), Y(1,2), ' H','Color',[0.5 0.5 0.5],'verticalAlignment', 'middle');
% text(Y(n1+n2+2,1), Y(n1+n2+2,2), ' S','Color','black','verticalAlignment', 'middle');
% 
% set(gcf,'color','w');
% hold off;
% export_fig xavierPcaSick.pdf



% %% PCA using all Patients from healthy mean
% clear all;
% load('xavier_splines_v1_N30_noloops','dPatients_noloops');
% splineData = constructEmptySplineData;
% splineData.N = 30; %no. control points, must be bigger than n+1
% splineData.Nt = 10 + 2; %Number of time control points
% splineData.Nphi = 6; %No. control points for diffeomorphisms
% splineData.nS = 3; %spacial degree
% splineData.nT = 2; %time degree
% splineData.nPhi = 3; %diffemorphism degree
% splineData.quadDegree = [8,4]; %Quadrature precission
% splineData.dSpace = 2;
% splineData.noInterpolS = 2 * splineData.N; % For composition
% splineData = constructKnots(splineData);
% [quadData, quadDataTensor] = setupQuadData(splineData);
% dSick = {dPatients_noloops{10:19}};
% dHealthy = {dPatients_noloops{1:9}};
% splineData.a = [1 1 1];
% close all;
% noCurves=19;
% noAllCurves=21;
% dSpace=2;
% n1 = 9;
% n2 = 10;
% N=splineData.N;
% load('InitialVelHealthy.mat');
% load('InitialVelMeanToMean.mat');
% load('CHealthyA1B1C1.mat');
% load('CSickA1B1C1.mat');
% V = zeros(splineData.N*splineData.dSpace,noCurves);
% for ii = 1:noCurves;
%     V(:,ii) = InitialVelHealthy{ii}(:);
% end
% G = metricMatrixH2(CHealthy,splineData,quadData);
% Gsqrt = sqrtm(G);
% Sigma = zeros([N*dSpace, N*dSpace]);
% for jj = 1:noCurves;
%     Sigma = Sigma + V(:,jj) * V(:,jj)';
% end
% 
% 
% Sigma = 1./(noCurves-1) * Gsqrt * Sigma * Gsqrt;
% 
% [U, Lambda] = eig(Sigma);
% Lambda = real(diag(Lambda));
% 
% % Calculate 2d coordinates (also of the means)
% for ii = 1:19;
%     V(:,ii) = InitialVelHealthy{ii}(:);
% end
% W = U(:,1:2);
% pts2d = zeros([2, 21]);
% Vproj = zeros([2, 21]);
% sz = size(V(:,1));
% VHealthy=zeros(sz);
% for jj = 2:n1+1
%     pts2d(:,jj) = W' * Gsqrt * V(:,jj-1);
%     Vproj(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
% end
% for jj = n1+2:n1+n2+1
%     pts2d(:,jj) = W' * Gsqrt * V(:,jj-1);
%     Vproj(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
% end
% VSick = -InitialVelMeanToMean{1}(:);
% 
% pts2d(:,1) = W' * Gsqrt * VHealthy;
% Vproj(:,1) = pts2d(:,1) ./ sqrt(Lambda(1:2));
% pts2d(:,noAllCurves) = W' * Gsqrt * VSick;
% Vproj(:,noAllCurves) = pts2d(:,noAllCurves) ./ sqrt(Lambda(1:2));
% 
% 
% %Plot
% hold on;
% axis equal;
% axis([-1 2 -1.7 3])
% grid on;
% box on;
% Y = Vproj';
% %plot(0,0,'o')
% %text(0, 0, ['mP'],'Color','black','verticalAlignment', 'middle');
% 
% plot(Y(2:n1+1,1),Y(2:n1+1,2), 'o','Color','blue');
% plot(Y(n1+2:n1+n2+1,1),Y(n1+2:n1+n2+1,2), 'x','Color','red');
% 
% plot(Y(1,1),Y(1,2), 'o','Color','black');
% plot(Y(n1+n2+2,1),Y(n1+n2+2,2), 'x','Color','red');
% for kk = 2:n1+1
%     text(Y(kk,1), Y(kk,2), [' ',num2str(kk-1)],'Color','blue','verticalAlignment', 'middle');
% end
% for kk = n1+2:n1+n2+1
%     text(Y(kk,1), Y(kk,2), [' ',num2str(kk-1)],'Color','red','verticalAlignment', 'middle');
% end
% daspect([.913  1  1])
% text(Y(1,1), Y(1,2), ' H','Color','black','verticalAlignment', 'middle');
% text(Y(n1+n2+2,1), Y(n1+n2+2,2), ' S','Color','red','verticalAlignment', 'middle');
% 
% set(gcf,'color','w');
% hold off;
% export_fig xavierPcaHealthyAllcurves.pdf


%% PCA using all Patients from sick mean
% clear all;
% load('xavier_splines_v1_N30_noloops','dPatients_noloops');
% splineData = constructEmptySplineData;
% splineData.N = 30; %no. control points, must be bigger than n+1
% splineData.Nt = 10 + 2; %Number of time control points
% splineData.Nphi = 6; %No. control points for diffeomorphisms
% splineData.nS = 3; %spacial degree
% splineData.nT = 2; %time degree
% splineData.nPhi = 3; %diffemorphism degree
% splineData.quadDegree = [8,4]; %Quadrature precission
% splineData.dSpace = 2;
% splineData.noInterpolS = 2 * splineData.N; % For composition
% splineData = constructKnots(splineData);
% [quadData, quadDataTensor] = setupQuadData(splineData);
% dSick = {dPatients_noloops{10:19}};
% dHealthy = {dPatients_noloops{1:9}};
% splineData.a = [1 1 1];
% noCurves=19;
% noAllCurves=21;
% dSpace=2;
% n1 = 9;
% n2 = 10;
% N=splineData.N;
% load('InitialVelSick.mat');
% load('InitialVelMeanToMean.mat');
% load('CHealthyA1B1C1.mat');
% load('CSickA1B1C1.mat');
% V = zeros(splineData.N*splineData.dSpace,noCurves);
% for ii = 1:noCurves;
%     V(:,ii) = InitialVelSick{ii}(:);
% end
% G = metricMatrixH2(CSick,splineData,quadData);
% Gsqrt = sqrtm(G);
% Sigma = zeros([N*dSpace, N*dSpace]);
% for jj = 1:noCurves;
%     Sigma = Sigma + V(:,jj) * V(:,jj)';
% end
% 
% 
% Sigma = 1./(noCurves-1) * Gsqrt * Sigma * Gsqrt;
% 
% [U, Lambda] = eig(Sigma);
% Lambda = real(diag(Lambda));
% 
% % Calculate 2d coordinates (also of the means)
% for ii = 1:19;
%     V(:,ii) = InitialVelSick{ii}(:);
% end
% W = U(:,1:2);
% pts2d = zeros([2, 21]);
% Vproj = zeros([2, 21]);
% sz = size(V(:,1));
% VSick=zeros(sz);
% for jj = 2:n1+1
%     pts2d(:,jj) = W' * Gsqrt * V(:,jj-1);
%     Vproj(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
% end
% for jj = n1+2:n1+n2+1
%     pts2d(:,jj) = W' * Gsqrt * V(:,jj-1);
%     Vproj(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
% end
% VHealthy = InitialVelMeanToMean{1}(:);
% 
% pts2d(:,1) = W' * Gsqrt * VHealthy;
% Vproj(:,1) = pts2d(:,1) ./ sqrt(Lambda(1:2));
% pts2d(:,noAllCurves) = W' * Gsqrt * VSick;
% Vproj(:,noAllCurves) = pts2d(:,noAllCurves) ./ sqrt(Lambda(1:2));
% 
% 
% %Plot
% figure(2)
% hold on;
% axis equal;
% axis([-2.7 0.5 -2.2 2.5])
% grid on;
% box on;
% Y = Vproj';
% %plot(0,0,'o')
% %text(0, 0, ['mP'],'Color','black','verticalAlignment', 'middle');
% 
% plot(Y(2:n1+1,1),Y(2:n1+1,2), 'o','Color','blue');
% plot(Y(n1+2:n1+n2+1,1),Y(n1+2:n1+n2+1,2), 'x','Color','red');
% 
% plot(Y(1,1),Y(1,2), 'o','Color','blue');
% plot(Y(n1+n2+2,1),Y(n1+n2+2,2), 'x','Color','black');
% for kk = 2:n1+1
%     text(Y(kk,1), Y(kk,2), [' ',num2str(kk-1)],'Color','blue','verticalAlignment', 'middle');
% end
% for kk = n1+2:n1+n2+1
%     text(Y(kk,1), Y(kk,2), [' ',num2str(kk-1)],'Color','red','verticalAlignment', 'middle');
% end
% 
% text(Y(1,1), Y(1,2), ' H','Color','blue','verticalAlignment', 'middle');
% text(Y(n1+n2+2,1), Y(n1+n2+2,2), ' S','Color','black','verticalAlignment', 'middle');
% pbaspect();
% daspect();
% set(gcf,'color','w');
% hold off;
% export_fig xavierPcaSickAllcurves.pdf
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
