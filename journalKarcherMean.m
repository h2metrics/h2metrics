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
splineData.a = [1 0 1e-4];
splineData.stepsT = 50;

options = struct( 'optDiff', true, ...
                  'optTra', true, ...
                  'optRot', true, ...
                  'optShift', true, ...
                  'tolFun', 1e-12, ...
                  'tolX', 1e-12, ...
                  'display', 'off', ... % 'off', 'iter-detailed'
                  'maxIter', 300, ...
                  'karcherTolGradNorm', 1e-3, ...
                  'karcherMaxIter', 30 );

d1 = loadDataSet('basic', splineData, '', ...
                 'curves', 'prop3', 'noise', 0.);
% d2 = loadDataSet('basic', splineData, '', ...
%                  'curves', 'prop3', 'noise', 0.1);

%Create some curves with noise
% noiseLevel = 0.10;
% d2 = d1 + noiseLevel*d1.*rand( size(d1) );
% d3 = d1 + noiseLevel*d1.*rand( size(d1) );
% d4 = d1 + noiseLevel*d1.*rand( size(d1) );
% d5 = d1 + noiseLevel*d1.*rand( size(d1) );
% d6 = d1 + noiseLevel*d1.*rand( size(d1) );
% d7 = d1 + noiseLevel*d1.*rand( size(d1) );
% d8 = d1 + noiseLevel*d1.*rand( size(d1) );
% d9 = d1 + noiseLevel*d1.*rand( size(d1) );
             
% [Eopt, dPathOpt,gOpt,info] = geodesicBvp(d1,d2,splineData,'options',options);             
             
% dListNoise15 = {d2,d3,d4,d5,d6,d7,d8,d9};

%% Compute, uncomment to compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remember to have this file!!!!
load([prefixDir, 'karcherMeanProp3.mat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dList = meanSplinesTest1.dList;

%
[dMeanProp, info] = karcherMeanManopt( dList, ...
    splineData, 'options', options, 'meanInit', d1 );

% Save
splineDataSave = splineData;
splineDataSave = rmfield(splineDataSave, 'quadData');
splineDataSave = rmfield(splineDataSave, 'quadDataTensor');
meanSplinesTest1 = struct('dList',{dList},'dMean',dMeanProp,'splineData',splineDataSave,'options',options,'info',info);
save('karcherMeanProp3.mat','meanSplinesTest1');

%% Plot and save some figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remember to have this file!!!!
load([prefixDir, 'karcherMeanProp3.mat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dMeanProp = meanSplinesTest1.dMean;
info = meanSplinesTest1.info;
% plotDir = 'journal/karcher_mean/';

%Data curves
figure(1)
for i = 1:length(dList)
clf
plotCurve( dList{i}, splineData,'-k' )
axis off
set(gcf,'color','w');
export_fig([plotDir,'propellerMeanErrord',num2str(i),'.pdf']);    
end

%Mean
clf
plotCurve( dMeanProp, splineData,'-k' )
axis off
set(gcf,'color','w');
export_fig([plotDir,'propellerMeanError.pdf']);

%Iterations in Karcher Mean algorithm
for j = [1,2,3,4,5,10,15] %Which iterations to plot
clf
plotCurve( info(j+1).dIter, splineData,'-k' )
axis off
set(gcf,'color','w');
export_fig([plotDir,'propellerMeanIter',num2str(j),'.pdf']);    
end

%% Big plot.
% plotDir = 'journal/final_plots/';
load([prefixDir, 'karcherMeanProp3.mat']);

dMean = meanSplinesTest1.dMean;

splineDataPlot = meanSplinesTest1.splineData;
splineDataPlot = setupQuadData(splineDataPlot);


dList = {};
for ii = 1:length(meanSplinesTest1.dList)
    dList{ii} = curveCenter(meanSplinesTest1.dList{ii}, splineDataPlot );
end

%% Start plotting 4x2 curves
handle = figure( 'Units', 'centimeters', 'Position', [0, 0, 12 4], ...
                 'Color', 'white' );
handle.Visible = 'off';
clf;
hold on;
axis equal;
axis off;

noPlotPtsS = 100;
plotPtsS = linspace(0, 2*pi, noPlotPtsS);

xmin = 0;
xmax = 0;
ymin = 0;
ymax = 0;

ax = [];

for jj = 2:-1:1
    for kk = 4:-1:1
        ax(4*(jj-1) + kk) = axes( 'Clipping', 'off', 'Position', ...
                                  [(kk-1)/6, (jj-1)/2, 1/6, 1/2] );
        c = deBoor( splineDataPlot.knotsS, splineDataPlot.nS,...
                dList{4*(jj-1)+kk}, plotPtsS, 1, ...
                    'periodic', true );
        plot(c(:,1), c(:,2), 'k-', 'Clipping', 'off');
        
        xl = xlim();
        xmin = min(xmin, xl(1));
        xmax = max(xmax, xl(2));
        xlim([xmin, xmax]);
        
        yl = ylim();
        ymin = min(ymin, yl(1));
        ymax = max(ymax, yl(2));
        ylim([ymin, ymax]);
        
        axis equal;
        axis off;
        axis tight;
    end
end
xl = xlim();
xlim([xl(1)-0.05*(xl(2)-xl(1)), xl(2)+0.05*(xl(2)-xl(1))]);
yl = ylim();
ylim([yl(1)-0.05*(yl(2)-yl(1)), yl(2)+0.05*(yl(2)-yl(1))]);

linkaxes(ax);

% Now plot the mean(s)
axMean = axes( 'Clipping', 'off', ...
               'Position', [4/6, 0, 2/6, 1]);
c = deBoor( splineDataPlot.knotsS, splineDataPlot.nS, dMean, plotPtsS, 1, ...
            'periodic', true );
plot(c(:,1), c(:,2), 'k-', 'Clipping', 'off');
axis equal;
axis off;
axis tight;

xl = xlim();
xlim([xl(1)-0.05*(xl(2)-xl(1)), xl(2)+0.05*(xl(2)-xl(1))]);
yl = ylim();
ylim([yl(1)-0.05*(yl(2)-yl(1)), yl(2)+0.05*(yl(2)-yl(1))]);

% Finish up plotting
hold off;

figname = [ plotDir, 'propellersMeanBig.pdf' ];
export_fig(figname);



             