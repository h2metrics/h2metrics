function plotData = setupPlotData( plotPointsS,plotPointsT,spData)
%UNTITLED5 
%
% TODO: Summary

%For plotting
% Input: d, t-values in [0,1] for plotting, s-values in [0,1] for plotting
% Output: Plot?

noPlotPointsS = length(plotPointsS);
noPlotPointsT = length(plotPointsT);

nS = spData.nS;
nT = spData.nT;
N = spData.N;
Nt = spData.Nt;
% quadDegree = spData.quadDegree;

knotsS =  spData.knotsS;
knotsT = spData.knotsT;
% innerKnotsS = spData.innerknotsS;
% innerKnotsT = spData.innerknotsT;


%Plotting, derivatives of order [0,1,2]x[0,1]
B_S_plot = spcol( knotsS, nS+1, brk2knt( plotPointsS, 3 ),'sparse');
B_T_plot = spcol( knotsT, nT+1, brk2knt( plotPointsT, 2 ),'sparse');
%Periodic B-splines
B_S_plot_per = [B_S_plot(:,1:nS) + B_S_plot(:,end-nS+1:end), B_S_plot(:,nS+1:end-nS)];

plotData = struct('points',[],'noPlotPointsS',noPlotPointsS,...
    'noPlotPointsT',noPlotPointsT,'B',[],'Bu',[],'Bt',[],'Buu',[],...
    'But',[],'Buut',[]);

plotData.points = {plotPointsS,plotPointsT};

for ii = N:-1:1
    for jj = Nt:-1:1;
     
    plotData.B(:,ii+(jj-1)*N) = reshape(B_S_plot_per(1:3:end,ii)*B_T_plot(1:2:end,jj)',...
        noPlotPointsS*noPlotPointsT ,1) ;
    plotData.Bu(:,ii+(jj-1)*N) = reshape(B_S_plot_per(2:3:end,ii)*B_T_plot(1:2:end,jj)',...
        noPlotPointsS*noPlotPointsT ,1)' ;
    plotData.Bt(:,ii+(jj-1)*N) = reshape(B_S_plot_per(1:3:end,ii)*B_T_plot(2:2:end,jj)',...
        noPlotPointsS*noPlotPointsT ,1)' ;
    plotData.Buu(:,ii+(jj-1)*N) = reshape(B_S_plot_per(3:3:end,ii)*B_T_plot(1:2:end,jj)',...
        noPlotPointsS*noPlotPointsT ,1)' ;
    plotData.But(:,ii+(jj-1)*N) = reshape(B_S_plot_per(2:3:end,ii)*B_T_plot(2:2:end,jj)',...
        noPlotPointsS*noPlotPointsT ,1)' ;
    plotData.Buut(:,ii+(jj-1)*N) = reshape(B_S_plot_per(3:3:end,ii)*B_T_plot(2:2:end,jj)',...
        noPlotPointsS*noPlotPointsT ,1)' ;
    
    end
end

end

