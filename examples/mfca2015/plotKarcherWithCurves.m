function  plotKarcherWithCurves(dKarcher,d, splineData)
% Produce a plot of a set of spline curves
n = length(d);
%% Plot results, H2
noPlotPointsS = 100;
plotPointsS = linspace(0,2*pi,noPlotPointsS);
%noPlotPointsT = 10;
plotPointsT = linspace(0,1,1);
plotData = setupPlotData(plotPointsS,plotPointsT,splineData);
figure;
hold all
ColorSet = varycolor(n);
set(gca, 'ColorOrder', ColorSet);
set(gcf,'color','w');
axis equal;
for i=1:n
    d_path =  linearPath(d{i},d{i},splineData);
    C = plotData.B*d_path;
    plot(C(1:plotData.noPlotPointsS,1),C(1:plotData.noPlotPointsS,2),':');   
end
d_path =  linearPath(dKarcher,dKarcher,splineData);
C = plotData.B*d_path;
plot(C(1:plotData.noPlotPointsS,1),C(1:plotData.noPlotPointsS,2),'k','LineWidth',2);
hold off
end

