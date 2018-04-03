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

%% Plot Karcher means
%A1B1C1
load('CSickA1B1C1');
plotKarcherWithCurves(CSick,dSick,splineData)
ylim(ylim-[.1,0])
export_fig xavierCSick.pdf

load('CHealthyA1B1C1');
plotKarcherWithCurves(CHealthy,dHealthy,splineData)
ylim(ylim-[.1,0])
export_fig xavierCHealthy.pdf

load('CPatientsA1B1C1');
plotKarcherWithCurves(CPatients,dPatients_noloops,splineData)
ylim(ylim-[.1,0])
export_fig xavierCPatients.pdf


%% Plot Geodesic between Karchermeans
%A1B1C1
load('dPathA1B1C1MeanToMean');
plotParticlePath(dPathA1B1C1MeanToMean,splineData,6,40);
export_fig xavierMean2Mean.pdf
plotPathRow(dPathA1B1C1MeanToMean,splineData,7);
export_fig xavierMean2MeanRow.pdf



%% Dendogramm Plot
load('DistA1B1C1Means.mat');
Dist = DistA1B1C1Means(2:20,2:20);
Y = squareform(Dist,'tovector');
Z = linkage(Y);
c = cluster(Z,'maxclust',3);
set(gcf,'color','w');
dendrogram(Z)

export_fig xavierDendogramm.pdf



%% Plot Mds of distance matrix
%A1B1C1withMeans
load('DistA1B1C1Means.mat');
tmpDist = DistA1B1C1Means;
[Y, e] = cmdscale(tmpDist);
close all
hold on;
n1=10;
n2=length(Y(:,1));
%deltax=0.02;
%deltay=-0.007;
plot(Y(1,1),-Y(1,2), 'o','Color','black');
plot(Y(2:n1,1),-Y(2:n1,2), 'o','Color','blue');
plot(Y(n1+1:n2-1,1),-Y(n1+1:n2-1,2), 'x','Color','red');
plot(Y(n2,1),-Y(n2,2), 'x','Color','black');
hold on;
axis equal;
grid on;
box on;
for kk = 2:n1
    text(Y(kk,1), -Y(kk,2), [' ',num2str(kk-1)],'Color','blue','verticalAlignment', 'middle');
end
for kk = n1+1:n2-1
    text(Y(kk,1), -Y(kk,2), [' ',num2str(kk-1)],'Color','red','verticalAlignment', 'middle');
end
text(Y(1,1), -Y(1,2), ' mH','Color','black','verticalAlignment', 'middle');
text(Y(n2,1), -Y(n2,2), ' mS','Color','black','verticalAlignment', 'middle');

set(gcf,'color','w');
axis([-0.6 0.6 -0.6 0.35]);
hold off;
export_fig xavierMdsMeans.pdf



%% Plot Mds of distance matrix without means 
%A1B1C1withMeans
load('DistA1B1C1Means.mat');
tmpDist = DistA1B1C1Means(2:20,2:20);
[Y, e] = cmdscale(tmpDist);
close all
hold on;
n1=9;
n2=10;
%deltax=0.02;
%deltay=-0.007;
plot(Y(1:n1,1),Y(1:n1,2), 'o','Color','blue');
plot(Y(n1+1:n1+n2,1),Y(n1+1:n1+n2,2), 'x','Color','red');
hold on;
axis equal;
grid on;
box on;
for kk = 1:n1
    text(Y(kk,1), Y(kk,2), [' ',num2str(kk)],'Color','blue','verticalAlignment', 'middle');
end
for kk = n1+1:n1+n2
    text(Y(kk,1), Y(kk,2), [' ',num2str(kk)],'Color','red','verticalAlignment', 'middle');
end
set(gcf,'color','w');
axis([-0.6 0.6 -0.3 0.65]);
hold off;
export_fig xavierMds.pdf





%% PCA Analysis of initial velocites to mean
% load('InitialVelPatientsA1B1C1.mat');
% load('CPatientsA1B1C1.mat');
% V = zeros(splineData.N*splineData.dSpace,length(InitialVelPatientsA1B1C1));
% for ii = 1:length(InitialVelPatientsA1B1C1);
%     V(:,ii) = InitialVelPatientsA1B1C1{ii}(:);
% end
% G = metricMatrixH2( CPatients ,splineData,quadData);
% Gsqrt = sqrtm(G);
% Vsqrt = Gsqrt*V;
% %Vpca = pca(Vsqrt');
% [Vpca2,Vpca2Score,latent]  = princomp(Vsqrt');
% Vprin = Vpca2'*Vsqrt;
% Vproj = Vprin(1:2,:);
% 
% %Plot
% plot(Vproj(1,1:9),Vproj(2,1:9),'o',Vproj(1,10:end),Vproj(2,10:end),'rx')
% axis equal
% axis off
% n1 = 9;
% n2 = 10;
% Y = Vproj';
% for kk = 1:n1
%     text(Y(kk,1), Y(kk,2), [' ',num2str(kk)],'Color','blue','verticalAlignment', 'middle');
% end
% for kk = n1+1:n1+n2
%     text(Y(kk,1), Y(kk,2), [' ',num2str(kk)],'Color','red','verticalAlignment', 'middle');
% end
% set(gcf,'color','w');
% export_fig xavierPca.pdf










%% Other constants in the metric

%A50B30C20
% load('CSickA50B30C20');
% plotKarcherWithCurves(CSick,dSick,splineData)
% ylim(ylim-[.1,0])
% export_fig CSickA50B30C20.pdf
% 
% load('CHealthyA50B30C20');
% plotKarcherWithCurves(CHealthy,dHealthy,splineData)
% ylim(ylim-[.1,0])
% export_fig CHealthyA50B30C20.pdf
% 
% %A20B30C50
% load('CSickA20B30C50');
% plotKarcherWithCurves(CSick,dSick,splineData)
% ylim(ylim-[.1,0])
% export_fig CSickA20B30C50.pdf
% 
% load('CHealthyA20B30C50');
% plotKarcherWithCurves(CHealthy,dHealthy,splineData)
% ylim(ylim-[.1,0])
% export_fig CHealthyA20B30C50.pdf
% 
% %A10B10C80
% load('CSickA10B10C80');
% plotKarcherWithCurves(CSick,dSick,splineData)
% ylim(ylim-[.1,0])
% export_fig CSickA10B10C80.pdf
% 
% load('CHealthyA10B10C80');
% plotKarcherWithCurves(CHealthy,dHealthy,splineData)
% ylim(ylim-[.1,0])
% export_fig CHealthyA10B10C80.pdf


%  load('dPathA80B10C10MeanToMean');
%  plotParticlePath(dPath,splineData,6,40);
%  export_fig Mean2MeanA80B10C10.pdf
%  plotPathRow(dPath,splineData,6);
%  export_fig Mean2MeanA80B10C10Row.pdf

























