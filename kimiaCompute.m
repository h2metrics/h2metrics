%% Load curves
load('kimiasplinesN40ns3');
extractnames = fieldnames(rmfield(kimiaControlPoints,{'N','nS'}));
for ii = 1:length(extractnames)
    dShapes{ii} = kimiaControlPoints.(extractnames{ii});
end
noCurves = length(dShapes);
splineData = constructEmptySplineData;
splineData.N = 40; %no. control points, must be bigger than n+1
splineData.Nt = 10 + 2; %Number of time control points
splineData.Nphi = 6; %No. control points for diffeomorphisms
splineData.nS = 3; %spacial degree
splineData.nT = 2; %time degree
splineData.nPhi = 3; %diffemorphism degree
splineData.quadDegree = [8,4]; %Quadrature precission
splineData.dSpace = 2;
splineData.noInterpolS = 2 * splineData.N; % For composition
splineData = constructKnots(splineData);
dSpace = 2;
[quadData, quadDataTensor] = setupQuadData(splineData);
splineData.a = [1 1 1];

%%
%CAT to COW (load curves and rotate)
dCurves = {kimiaControlPoints.cat1,kimiaControlPoints.cow1};
%Rotate and scale curves
A=zeros(2);
A(1,2)=1;
A(2,1)=-1;
noCurves = length(dCurves);
for kk = noCurves:-1:1
    cQuad = quadData.B_S * dCurves{kk};
    cQuad_u = quadData.Bu_S * dCurves{kk};
    cSpeed = sum( cQuad_u.^2 , 2).^(1/2); 
    cLength = sum(cSpeed .* quadData.quadWeightsS);
    dCurves{kk} = A*(dCurves{kk}');
    dCurves{kk} = dCurves{kk}';
    dCurves{kk} = dCurves{kk} / cLength;
end
[A B C] = determineConstants(dCurves,splineData,quadData,quadDataTensor);
splineData.a = [A B C];
[E_geo, dPath1] = geodesicBvpAmpl(dCurves{1},dCurves{2},splineData,quadData,quadDataTensor);
save('Cat2CowA1B1C1','dPath1');
splineData.a = [10*A 10*B 80*C];
[E_geo, dPath2] = geodesicBvpAmpl(dCurves{1},dCurves{2},splineData,quadData,quadDataTensor);
save('Cat2CowA10B10C80','dPath2');
splineData.a = [1*A 1*B 98*C];
[E_geo, dPath3] = geodesicBvpAmpl(dCurves{1},dCurves{2},splineData,quadData,quadDataTensor);
save('Cat2CowA1B1C98','dPath3');


%% Plot
noPlotPointsS = 100;
noPlotPointsT = 6;
plotPointsS = linspace(0,2*pi,noPlotPointsS);
plotPointsT = linspace(0,1,noPlotPointsT);
plotData = setupPlotData(plotPointsS,plotPointsT, splineData);
c1 = plotData.B*dPath1;
c2 = plotData.B*dPath2;
c3 = plotData.B*dPath3;

x_size1 = (max(max(c1(:,1)))-min(min(c1(:,1))));
x_max1= max(max(c1(:,1)))+x_size1/20;
x_min1= min(min(c1(:,1)))-x_size1/20;


x_size2 = (max(max(c2(:,1)))-min(min(c2(:,1))));
x_max2= max(max(c2(:,1)))+x_size2/20;
x_min2= min(min(c2(:,1)))-x_size2/20;


x_size3 = (max(max(c3(:,1)))-min(min(c3(:,1))));
x_max3= max(max(c3(:,1)))+x_size3/20;
x_min3= min(min(c3(:,1)))-x_size3/20;

x_max = max([x_max1 x_max2 x_max3]);
x_min = max([x_min1 x_min2 x_min3]);


figure;
hold all;
set(gcf,'color','w');
axis off;
axis equal;
for ii = 1:plotData.noPlotPointsT
    plot(c1( 1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 1)+1.1*(ii-1)*(x_max-x_min),...
        c1(1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 2),'color','black');
end
export_fig test1.pdf

hold off;

figure;
hold all;
set(gcf,'color','w');
axis off;
axis equal;
for ii = 1:plotData.noPlotPointsT
    plot(c2( 1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 1)+1.1*(ii-1)*(x_max-x_min),...
        c2(1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 2),'color','black');
end
export_fig test2.pdf

hold off;


figure;
hold all;
set(gcf,'color','w');
axis off;
axis equal;
for ii = 1:plotData.noPlotPointsT
    plot(c3( 1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 1)+1.1*(ii-1)*(x_max-x_min),...
        c3(1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 2),'color','black');
end
export_fig test3.pdf

hold off;









