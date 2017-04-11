function   plotPathRow(dPath,splineData,noPlotPointsT)
noPlotPointsS = 100;
plotPointsS = linspace(0,2*pi,noPlotPointsS);
plotPointsT = linspace(0,1,noPlotPointsT);
plotData = setupPlotData(plotPointsS,plotPointsT, splineData);
c = plotData.B*dPath;

x_size = (max(max(c(:,1)))-min(min(c(:,1))));
x_max= max(max(c(:,1)))+x_size/20;
x_min= min(min(c(:,1)))-x_size/20;

figure;
hold all;
set(gcf,'color','w');
ColorSet = varycolor((noPlotPointsT));
set(gca, 'ColorOrder', ColorSet);
axis off;
axis equal;
for ii = 1:plotData.noPlotPointsT
    plot(c( 1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 1)+1.1*(ii-1)*(x_max-x_min),...
        c(1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 2) );
end
%export_fig test.pdf

hold off;










end   