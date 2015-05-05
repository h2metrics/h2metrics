function  plotPath( dPath, splineData, noPlotPointsT)
% Produce a plot of a spline path

%% Plot results, H2
noPlotPointsS = 100;
plotPointsS = linspace(0,2*pi,noPlotPointsS);
%noPlotPointsT = 10;
plotPointsT = linspace(0,1,noPlotPointsT);
plotData = setupPlotData( plotPointsS,plotPointsT, splineData);

% %% For plotting
% C =quadData.B*d_energyH2;
% To plot the basis function
%surf( quadData.points{2}', quadData.points{1}', reshape(quadData.B(:,3), 32,25) )

figure(2)
%dPath_optimal = [d0;d_optimal;d1];
%d_optimal_full = d_linear;

C_optimal = plotData.B*dPath;

clf
hold on
for ii = 1:plotData.noPlotPointsT
    plot( C_optimal( 1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 1),...
        C_optimal( 1 + (ii-1)*plotData.noPlotPointsS: (ii-1)*plotData.noPlotPointsS + plotData.noPlotPointsS, 2) );
    axis equal
end
hold off

%% Figure 3, snapshots
noSnapshots = 9;
plotSnapShots = linspace(0,1,noSnapshots);
plotDataSH = setupPlotData( plotPointsS,plotSnapShots, splineData);
C_SH = plotDataSH.B*dPath;
figure(3)
for ii = 1:noSnapshots
    subplot( 3,3, ii)
    plot( C_SH( 1 + (ii-1)*plotDataSH.noPlotPointsS: (ii-1)*plotDataSH.noPlotPointsS + plotDataSH.noPlotPointsS, 1),...
        C_SH( 1 + (ii-1)*plotDataSH.noPlotPointsS: (ii-1)*plotDataSH.noPlotPointsS + plotDataSH.noPlotPointsS, 2) );
    axis equal
end


end

