function   plotPathSingleThickThin(dPath,splineData,noPlotPointsT)
noPlotPointsS = 1000;
plotPointsS = linspace(0,2*pi,noPlotPointsS);
plotPointsT = linspace(0,1,noPlotPointsT);
plotData = setupPlotData(plotPointsS,plotPointsT, splineData);
c = plotData.B*dPath;

x_size = (max(max(c(:,1)))-min(min(c(:,1))));
y_size = (max(max(c(:,2)))-min(min(c(:,2))));
x_max= max(max(c(:,1)))+x_size/20;
x_min= min(min(c(:,1)))-x_size/20;
y_max= max(max(c(:,2)))+y_size/20;
y_min= min(min(c(:,2)))-y_size/20;

figure;
hold all;
set(gcf,'color','w');
%ColorSet = varycolor((noPlotPointsT));
%set(gca, 'ColorOrder', ColorSet);
axis off;
%axis([x_min x_max y_min y_max]);
 plot( c( 1: plotData.noPlotPointsS, 1),c( 1:plotData.noPlotPointsS, 2),'k','LineWidth',2);
 plot( c( 1 + (plotData.noPlotPointsT-1)*plotData.noPlotPointsS: (plotData.noPlotPointsT-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 1),...
        c( 1 + (plotData.noPlotPointsT-1)*plotData.noPlotPointsS: (plotData.noPlotPointsT-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 2),'k','LineWidth',2);  
for ii = 2:plotData.noPlotPointsT-1
    plot( c( 1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 1),...
        c( 1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 2),'k-');
    axis equal
end
hold off;



end   