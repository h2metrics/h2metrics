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
export_fig CSickA1B1C1.pdf

load('CHealthyA1B1C1');
plotKarcherWithCurves(CHealthy,dHealthy,splineData)
export_fig CHealthyA1B1C1.pdf

%A80B10C10
load('CSickA80B10C10');
plotKarcherWithCurves(CSick,dSick,splineData)
export_fig CSickA80B10C10.pdf

load('CHealthyA80B10C10');
plotKarcherWithCurves(CHealthy,dHealthy,splineData)
export_fig CHealthyA80B10C10.pdf


%% Plot Geodesic between Karchermeans
%A1B1C1

load('dPathA1B1C1MeanToMean');
plotParticlePath(dPath,splineData,6,40);
export_fig Mean2MeanA1B1C1.pdf
plotPathRow(dPath,splineData,6);
export_fig Mean2MeanA1B1C1Row.pdf


load('dPathA80B10C10MeanToMean');
plotParticlePath(dPath,splineData,6,40);
export_fig Mean2MeanA80B10C10.pdf
plotPathRow(dPath,splineData,6);
export_fig Mean2MeanA80B10C10Row.pdf





%% Plot Mds of distance matrix
load('DistA1B1C1Means.mat');
tmpDist = DistA1B1C1Means;

[Y, e] = cmdscale(tmpDist);

%% Multidimensional scaling - 2d plot, length
%A1B1C1withMeans
load('DistA1B1C1Means.mat');
tmpDist = DistA1B1C1Means;
[Y, e] = cmdscale(tmpDist);

close all
hold on;
n1=10;
n2=length(Y(:,1));
deltax=0.02;
deltay=-0.007;
plot(Y(1:n1,1), Y(1:n1,2), 'o','Color','blue');
plot(Y(n1+1:n2,1), Y(n1+1:n2,2), 'x','Color','red');
hold on;
for kk = 2:n1
    text(Y(kk,1)+deltax, Y(kk,2)+deltay, num2str(kk-1),'Color','blue');
end
for kk = n1+1:n2-1
    text(Y(kk,1)+deltax, Y(kk,2)+deltay, num2str(kk),'Color','red');
end
text(Y(1,1)+deltax, Y(1,2)+deltay, 'd_H','Color','blue');
text(Y(n2,1)+deltax, Y(n2,2)+deltay, 'd_S','Color','red');

set(gcf,'color','w');
hold off;
export_fig MdsA1B1C1Means.pdf

%A1B1C1withoutMeans
load('DistA1B1C1Means.mat');
tmpDist = DistA1B1C1Means;
[Y, e] = cmdscale(tmpDist);

close all
hold on;
n1=10;
n2=length(Y(:,1));
deltax=0.02;
deltay=-0.007;
plot(Y(2:n1,1), Y(2:n1,2), 'o','Color','blue');
plot(Y(n1+1:n2-1,1), Y(n1+1:n2-1,2), 'x','Color','red');
hold on;
for kk = 2:n1
    text(Y(kk,1)+deltax, Y(kk,2)+deltay, num2str(kk-1),'Color','blue');
end
for kk = n1+1:n2-1
    text(Y(kk,1)+deltax, Y(kk,2)+deltay, num2str(kk),'Color','red');
end
set(gcf,'color','w');
hold off;
export_fig MdsA1B1C1.pdf







