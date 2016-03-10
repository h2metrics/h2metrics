%% Plot and save figure for basic_figure
% load('journal/karcher_mean/karcherMeanProp3.mat');
% plotDir = 'journal/final_plots/';

splineData = constructEmptySplineData;
splineData.N = 70;
splineData.nS = 3;
splineData.Nt = 20;
splineData.nT = 2;
splineData.Nphi = 20;
splineData.nPhi = 3;
splineData.quadDegree = [6, 4];
splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);

dCircle = loadDataSet('basic', splineData, '', ...
                 'curves', 'circle', 'noise', 0.);             
dProp3Noise = loadDataSet('basic', splineData, '', ...
                 'curves', 'prop3', 'noise', 0.1);
dProp4Noise = loadDataSet('basic', splineData, '', ...
                 'curves', 'prop4', 'noise', 0.1);

             
dWrap = loadDataSet('basic', splineData, '', ...
                 'curves', 'wrap', 'noise', 0.);
dProp3 = loadDataSet('basic', splineData, '', ...
                 'curves', 'prop3', 'noise', 0.);
dProp4 = loadDataSet('basic', splineData, '', ...
                 'curves', 'prop4', 'noise', 0.);

dListTilak = loadDataSet( 'corpus_callosum_tilak', splineDataInit, dataDir, ...
                     'class', 'normal', 'ind', [1, 2], 'constSpeed');

lengthTilak = [ curveLength(dListTilak{1}, splineData),...
    curveLength(dListTilak{1}, splineData)];
dListTilakNorm = {dListTilak{1}/lengthTilak(1)*2*pi,...
    dListTilak{2}/lengthTilak(2)*2*pi}; 
                 
dList = {dCircle,dWrap,dProp3,dProp3Noise,dProp4,dProp4Noise,dListTilakNorm{1},dListTilakNorm{2}};

for ii = 1:length(dList)
    dList{ii} = curveCenter(dList{ii}, splineData );
end

%% Generate seperate pdf's

noPlotPtsS = 200;
plotPtsS = linspace(0, 2*pi, noPlotPtsS);

for ii = 1:length(dList)-2
    c = deBoor( splineDataPlot.knotsS, splineDataPlot.nS,...
                dList{ii}, plotPtsS, 1, 'periodic', true );
            
        plot(c(:,1), c(:,2), 'k-', 'Clipping', 'off','LineWidth',2);
        axis equal;
        axis off;
        axis tight;
set(gcf,'Color', 'white');
figname = [ plotDir, 'basicFiguresC',num2str(ii),'.pdf' ];
export_fig(figname);
end
%% Rotate corpus callosum, i = 4,8
for ii = 7:8
    c = deBoor( splineDataPlot.knotsS, splineDataPlot.nS,...
                dList{ii}, plotPtsS, 1, 'periodic', true );
        rotMatrix = [cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4) ];
        c = c*rotMatrix;
        plot(c(:,1), c(:,2), 'k-', 'Clipping', 'off','LineWidth',2);
        axis equal;
        axis off;
        axis tight;
set(gcf,'Color', 'white');
figname = [ plotDir, 'basicFiguresC',num2str(ii),'.pdf' ];
export_fig(figname);
end




