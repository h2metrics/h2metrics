%% Section 1: Setup 
dataDir = '..\data'; %Subdirectory for data


%% Load curves
splineData = constructEmptySplineData;
splineData.N = 60; %no. control points, must be bigger than n+1
splineData.Nt = 20 + 2; %Number of time control points
splineData.Nphi = 25; %No. control points for diffeomorphisms
splineData.nS = 3; %spacial degree
splineData.nT = 2; %time degree
splineData.nPhi = 3; %diffemorphism degree
splineData.quadDegree = [8,4]; %Quadrature precission
splineData.dSpace = 2;
splineData.noInterpolS = 2 * splineData.N; % For composition
splineData = constructKnots(splineData);
dSpace = 2;
splineData = setupQuadData(splineData);
%n = length(dShapes);
splineData.a = [1 0 1];

options = struct( 'optDiff', true, ...
                  'optTra', false, ...
                  'optRot', false, ...
                  'optShift', true, ...
                  'tolFun', 1e-12, ...
                  'tolX', 1e-12, ...
                  'display', 'iter-detailed', ... % 'off', 'iter-detailed'
                  'maxIter', 1000 ,...
                  'rigidGlobalRot', false, ...
                  'rigidUseComp', true);                         
         

          
d1 = loadDataSet('basic', splineData, '', ...
                 'curves', 'circle', 'noise', 0.);
d2 = loadDataSet('basic', splineData, '', ...
                 'curves', 'wrap', 'noise', 0.);


             

%% Optimize

[optE, optP1,optGamma] = geodesicBvp(d1,d2,splineData,'options', options);


%% Plot
noMarkers= 20;
dPath = optP1;
noPlotPointsS = 1000;
plotPointsS = linspace(0,2*pi,noPlotPointsS);
plotMarkerPoints =  linspace(0,2*pi,noMarkers);
plotPointsT = linspace(0,1,5);
plotData = setupPlotData(plotPointsS,plotPointsT, splineData);
plotData2 = setupPlotData(plotMarkerPoints,plotPointsT, splineData);

c = plotData.B*dPath;
quadData = splineData.quadData;
quadDataTensor = splineData.quadDataTensor;

%% Evaluate path at quadrature sites
Cu = quadDataTensor.Bu * dPath;
Ct = quadDataTensor.Bt * dPath;
Cut = quadDataTensor.But * dPath;
Cuu = quadDataTensor.Buu * dPath;
Cuut = quadDataTensor.Buut * dPath;





c2 = plotData2.B*dPath;

x_size = (max(max(c(:,1)))-min(min(c(:,1))));
x_max= max(max(c(:,1)))+x_size/20;
x_min= min(min(c(:,1)))-x_size/20;

figure;
hold all;
set(gcf,'color','w');
axis off;
axis equal; 
for ii = 1:plotData.noPlotPointsT
    plot(c( 1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 1)+1.1*(ii-1)*(x_max-x_min),...
        c(1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 2),'k');
    plot(c2( 1 + (ii-1)*plotData2.noPlotPointsS: (ii-1)*plotData2.noPlotPointsS + plotData2.noPlotPointsS, 1)+1.1*(ii-1)*(x_max-x_min),...
       c2(1 + (ii-1)*plotData2.noPlotPointsS: (ii-1)*plotData2.noPlotPointsS + plotData2.noPlotPointsS, 2),'ko','MarkerSize',4);
end
export_fig([plotDir, 'Forward.pdf'])
close
%plotDiff2(optGamma.phi, splineData)
%export_fig([plotDir,'/PhiForward.pdf'])
close

%% Optimize Backwards
[optEbw, optP1,optGamma] = geodesicBvp(d2,d1,splineData,'options', options);

%% Plot Backwards
dPath = optP1; 
noPlotPointsT = 5;
noPlotPointsS = 1000;
noPlotPointsS2 = 20;
plotPointsS = linspace(0,2*pi,noPlotPointsS);
plotPointsS2 = linspace(0,2*pi,noPlotPointsS2);
plotPointsT = linspace(0,1,noPlotPointsT);
plotData = setupPlotData(plotPointsS,plotPointsT, splineData);
plotData2 = setupPlotData(plotPointsS2,plotPointsT, splineData);
c = plotData.B*dPath;
c2 = plotData2.B*dPath;

x_size = (max(max(c(:,1)))-min(min(c(:,1))));
x_max= max(max(c(:,1)))+x_size/20;
x_min= min(min(c(:,1)))-x_size/20;

figure;
hold all;
set(gcf,'color','w');
axis off;
axis equal;
T = plotData.noPlotPointsT;
for ii = 1:T
    plot(c( 1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 1)+1.1*(T-ii)*(x_max-x_min),...
        c(1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 2),'k');
     plot(c2( 1 + (ii-1)*plotData2.noPlotPointsS: (ii-1)*plotData2.noPlotPointsS + plotData2.noPlotPointsS, 1)+1.1*(T-ii)*(x_max-x_min),...
        c2(1 + (ii-1)*plotData2.noPlotPointsS: (ii-1)*plotData2.noPlotPointsS + plotData2.noPlotPointsS, 2),'ok','MarkerSize',4);
end
hold off;
export_fig([plotDir,'Backward.pdf'])
close
%plotDiff2Inv(optGamma.phi, splineData)
%export_fig([plotDir,'/PhiBackward.pdf'])
close





    




