%% Script for convergence analysis
% HOW TO RUN: 
%Section 1: Setup a test case
%Section 2: Run experiments defined in section 1
%Section 3: Compute post processing results, H1H2 distances etc.
%Section 4: Plot results from section 3

%% Section 1: Setup Circle to wrap
% saveDir = 'journal/ConvergencePropellers/Data';
% dataDir = 'C:\Users\jakmo\Dropbox\H2_numerics\code\data'; %Subdirectory for data
% plotDir = 'journal/final_plots';
% dataName = 'ConvergencePropellers'; %name for the result files

%Define cases to calculate
NList = 10:10:110;
NtList = 10:5:90;
nSList = 3:5;
% nTList = 1:4;
nTList = 1:2;

%Construct a initial setup
splineDataInit = constructEmptySplineData;
splineDataInit.N = 70;
splineDataInit.nS = 3;
splineDataInit.Nt = 20;
splineDataInit.nT = 2;
splineDataInit.quadDegree = [6, 4];
splineDataInit = constructKnots(splineDataInit);
splineDataInit = setupQuadData(splineDataInit);
splineDataInit.a = [1 0 1e-4];

% aH1H2 = [splineDataInit.a;splineDataInit.a];
aH1H2 = [1,0,0;0,0,0];

%Load curve data
% d1Init = loadDataSet('basic', splineDataInit, '', ...
%                  'curves', 'circle');
% d2Init = loadDataSet('basic', splineDataInit, '', ...
%                  'curves', 'wrap');
d1 = loadDataSet('basic', splineDataInit, '', ...
                 'curves', 'prop3', 'noise', 0.);
d2 = loadDataSet('basic', splineDataInit, '', ...
                 'curves', 'prop4', 'noise', 0.);
             
             
%Setup optimization options
options = struct( 'optDiff', false, ...
                  'optTra', false, ...
                  'optRot', false, ...
                  'optShift', false, ...
                  'tolFun', 1e-12, ...
                  'tolX', 1e-12, ...
                  'display', 'iter-detailed', ... % 'off', 'iter-detailed'
                  'maxIter', 300 ,...
                  'rigidGlobalRot', false, ...
                  'rigidUseComp', true);
% Run optimization to find initial guess
[optEInit, optPathInit, optGaInit, infoInit] = geodesicBvp(d1, d2,...
    splineDataInit, ...
    'options', options);
% 
% [G comp] = pathRiemH2Energy(optPathInit,splineDataInit);
% comp(3)/G
% 
% figure(1)
% clf
% plotPath(optPathInit,splineDataInit)


%% Section 2: Compute experiments, takes ~6 minutes
numberOfExperiments = min( length(NList), length(NtList) );
energyDiagPlot = zeros(numberOfExperiments,length(nSList),length(nTList));
dofDiagPlot1 = zeros(numberOfExperiments,length(nSList),length(nTList));
pathDiag = cell(numberOfExperiments,length(nSList),length(nTList));

tic
for ii = 1:numberOfExperiments
    for jj = 1:length(nSList)
        for kk = 1:length(nTList)
            
            splineData = constructEmptySplineData;
            splineData.quadDegree = [6, 4]; %Quadrature precission
            splineData.dSpace = 2;
            splineData.a = [1 0 1e-4];
            
            splineData.N = NList(ii); %no. control points, must be bigger than n+1
            splineData.nS = nSList(jj); %spacial degree
            splineData.Nt = NtList(ii); %Number of time control points
            splineData.nT = nTList(kk); %time degree
            
            splineData = constructKnots(splineData);
            splineData = setupQuadData(splineData);
            
            d1 = loadDataSet('basic', splineData, '', ...
                'curves', 'prop3', 'noise', 0.);
            d2 = loadDataSet('basic', splineData, '', ...
                'curves', 'prop4', 'noise', 0.);
            
            
            %Initial guess, sample first solution
            dInit = pathSpline2Spline(optPathInit,splineDataInit,splineData);
            
            disp(['Case: N = ',num2str(NList(ii)),' ,Nt = ', num2str(NtList(ii)),...
                ' ,nS = ',num2str(nSList(jj)),' ,nT = ',num2str(nTList(kk))]);
            
            [optE, optPath, optGa, info] = geodesicBvp(d1, d2, splineData, ...
                'options', options); %'initpath',dInit);
            
%             savePath(dataName,optPath,splineData,...
%                 'workdir',saveDir,'energy',optE);
            
            %Save results
            energyDiagPlot(ii,jj,kk) = optE;
            dofDiagPlot1(ii,jj,kk) = splineData.N*splineData.Nt;
            pathDiag{ii,jj,kk} = optPath;
            
        end       
    end
end
toc
%Relative energy
energyRelDiffDiag = abs(diff(energyDiagPlot))./energyDiagPlot(2:end,:,:);

%% Plot relative energy difference of diagonals (for each nS and nT)
figure(1)
% for jj = 1:length(nSList)
%     for kk = 1:length(nTList)
%         loglog( dofDiagPlot1(2:end,jj,kk), energyRelDiffDiag(:,jj,kk),...
%             'Displayname',['nS=',num2str(nSList(jj)),' nT=',num2str(nTList(kk))] );
%         hold on
%     end
% end
loglog( dofDiagPlot1(1:end-1,1,1), energyRelDiffDiag(:,1,1),'-*k','LineWidth',1,...
            'Displayname',['$n_\theta$ = ',num2str(nSList(1)),', $n_t$ = ',num2str(nTList(1))] );
hold on
loglog( dofDiagPlot1(1:end-1,2,1), energyRelDiffDiag(:,2,1),'--*k','LineWidth',1,...
            'Displayname',['$n_\theta$ = ',num2str(nSList(2)),', $n_t$ = ',num2str(nTList(1))] );
loglog( dofDiagPlot1(1:end-1,3,1), energyRelDiffDiag(:,3,1),':*k','LineWidth',1,...
            'Displayname',['$n_\theta$ = ',num2str(nSList(3)),', $n_t$ = ',num2str(nTList(1))] );
loglog( dofDiagPlot1(1:end-1,1,2), energyRelDiffDiag(:,1,2),'-ok','LineWidth',1,...
            'Displayname',['$n_\theta$ = ',num2str(nSList(1)),', $n_t$ = ',num2str(nTList(2))] );
loglog( dofDiagPlot1(1:end-1,2,2), energyRelDiffDiag(:,2,2),'--ok','LineWidth',1,...
            'Displayname',['$n_\theta$ = ',num2str(nSList(2)),', $n_t$ = ',num2str(nTList(2))] );
loglog( dofDiagPlot1(1:end-1,3,2), energyRelDiffDiag(:,3,2),':ok','LineWidth',1,...
            'Displayname',['$n_\theta$ = ',num2str(nSList(3)),', $n_t$ = ',num2str(nTList(2))] );
hold off
% title('Relative energy difference on diagonals');
% ylabel('$|E_i - E_{i-1}|/E_i$', 'Interpreter','latex')
% ylabel('|E_i - E_{i-1}|/E_i')
% xlabel('$N_\theta \cdot N_t$', 'Interpreter','latex')
hleg = legend('show','location','southwest')
set(hleg,'Interpreter','latex','FontSize',12)
fig = gcf;
fig.Color = 'white';
export_fig([plotDir,'/RelEnergyDiagonalPropellers.pdf'])

%% Load "Diagonal Paths" and compute the relative H1H2 distance
dofDiagPlot = zeros(min(length(NList),length(NtList)),length(nSList),length(nTList));
pathH1H2DistDiag = zeros(min(length(NList),length(NtList))-1,length(nSList),length(nTList));
pathH1H2NormDiag = zeros(min(length(NList),length(NtList))-1,length(nSList),length(nTList));

pathH1H2DistDiagComp = zeros(6,min(length(NList),length(NtList))-1,length(nSList),length(nTList));
pathH1H2NormDiagComp = zeros(6,min(length(NList),length(NtList))-1,length(nSList),length(nTList));

%Initialize splineData
splineDataPrev = splineDataInit;
splineDataNext = splineDataInit;

for jj = 1:length(nSList)
    for kk = 1:length(nTList) 
        splineDataPrev.N = NList(1);
        splineDataPrev.Nt = NtList(1);
        splineDataPrev.nS = nSList(jj);
        splineDataPrev.nT = nTList(kk);
        splineDataPrev = constructKnots(splineDataPrev);
        [splineDataPrev] = setupQuadData(splineDataPrev);
        dPathPrev =  pathDiag{1,jj,kk};
        
        for ii = 2: min( length(NList), length(NtList) )
            splineDataNext.N = NList(ii);
            splineDataNext.Nt = NtList(ii);
            splineDataNext.nS = nSList(jj);
            splineDataNext.nT = nTList(kk);
            splineDataNext = constructKnots(splineDataNext);
            [splineDataNext] = setupQuadData(splineDataNext);
            dPathNext =  pathDiag{ii,jj,kk};
            
            
            disp(['Case: N = ',num2str(NList(ii)),' ,Nt = ', num2str(NtList(ii)),...
                     ' ,nS = ',num2str(nSList(jj)),' ,nT = ',num2str(nTList(kk))]);
            
            %Save information
%             pathH1H2DistDiag(ii-1,jj,kk) = pathFlatH1H2dist(dPathPrev,splineDataPrev,...
%                 dPathNext,splineDataNext,'a', aH1H2);
%             pathH1H2NormDiag(ii-1,jj,kk) = pathFlatH1H2Norm(dPathNext,splineDataNext,...
%                 'a', aH1H2);
%             
            %Calculate L2L2 norm difference
            [pathH1H2DistDiag(ii-1,jj,kk),pathH1H2DistDiagComp(:,ii-1,jj,kk)] = ...
                pathFlatH1H2Dist(dPathPrev,splineDataPrev,...
                dPathNext,splineDataNext,'a',aH1H2);
            [pathH1H2NormDiag(ii-1,jj,kk),pathH1H2NormDiagComp(:,ii-1,jj,kk)] = ...
                pathFlatH1H2Norm(dPathNext,splineDataNext,'a',aH1H2);

            dofDiagPlot(ii,jj,kk) = splineDataPrev.N*splineDataPrev.Nt;
            
            %Move to next spline
            dPathPrev = dPathNext;
            splineDataPrev = splineDataNext;
        end
    end
end
%Relative energy
pathH1H2RelDistDiag = pathH1H2DistDiag./pathH1H2NormDiag;

%% Plot diagonal H1H2 distances (for each nS and nT)
figure(2)
% for jj = 1:length(nSList)
%     for kk = 1:length(nTList)
%         loglog( dofDiagPlot(2:end,jj,kk), pathH1H2DistDiag(:,jj,kk),...
%             'Displayname',['nS=',num2str(nSList(jj)),' nT=',num2str(nTList(kk))] );
%         hold on
%     end
% end
% hold off
loglog( dofDiagPlot(2:end,1,1), pathH1H2DistDiag(:,1,1),'-xk','LineWidth',1,...
            'Displayname',['nS=',num2str(nSList(1)),' nT=',num2str(nTList(1))] );
hold on
loglog( dofDiagPlot(2:end,2,1), pathH1H2DistDiag(:,2,1),'--xk','LineWidth',1,...
            'Displayname',['nS=',num2str(nSList(2)),' nT=',num2str(nTList(1))] );
loglog( dofDiagPlot(2:end,3,1), pathH1H2DistDiag(:,3,1),':xk','LineWidth',1,...
            'Displayname',['nS=',num2str(nSList(3)),' nT=',num2str(nTList(1))] );
loglog( dofDiagPlot(2:end,1,2), pathH1H2DistDiag(:,1,2),'-ok','LineWidth',1,...
            'Displayname',['nS=',num2str(nSList(1)),' nT=',num2str(nTList(2))] );
loglog( dofDiagPlot(2:end,2,2), pathH1H2DistDiag(:,2,2),'--ok','LineWidth',1,...
            'Displayname',['nS=',num2str(nSList(2)),' nT=',num2str(nTList(2))] );
loglog( dofDiagPlot(2:end,3,2), pathH1H2DistDiag(:,3,2),':ok','LineWidth',1,...
            'Displayname',['nS=',num2str(nSList(3)),' nT=',num2str(nTList(2))] );
hold off
% set(gca,'Xlim',[0 10000],...
%     'XTick',[10^2 5*1e2 1e3 5*1e3 1e4] );
%         'YLim',[0,max(max(max( pathH1H2DistDiag )))],...
        
% title('Absolute H1H2 distance');
% ylabel('$||d_i - d_{i-1}||_{H1H2}$', 'Interpreter','latex')
% xlabel('DOF = N*Nt')
legend('show','location','northeast')
set(gcf, 'Color','white');
% export_fig ConvergenceCircle2Wrap/pdf/AbsH1H2DistDiagonal.pdf    
export_fig([plotDir,'/AbsL2L2DistDiagonalPropellers.pdf'])

%Relative distance
figure(3)
% for jj = 1:length(nSList)
%     for kk = 1:length(nTList)
%         loglog( dofDiagPlot(2:end,jj,kk), pathH1H2RelDistDiag(:,jj,kk),...
%             'Displayname',['nS=',num2str(nSList(jj)),' nT=',num2str(nTList(kk))] );
% %         legend(['nS=',num2str(nSList(jj)),'nT=',num2str(nTList(kk))])
%         hold on
% %         abs(diff(diag(energyArray(:,:,jj,kk,1))))./diag(energyArray(2:end,2:end,jj,kk,1))
%     end
% end
% hold off

hold off
loglog( dofDiagPlot(2:end,1,1), pathH1H2RelDistDiag(:,1,1),'-xk','LineWidth',1,...
            'Displayname',['$n_\theta$ = ',num2str(nSList(1)),', $n_t$ = ',num2str(nTList(1))] );
hold on
loglog( dofDiagPlot(2:end,2,1), pathH1H2RelDistDiag(:,2,1),'--xk','LineWidth',1,...
            'Displayname',['$n_\theta$ = ',num2str(nSList(2)),', $n_t$ = ',num2str(nTList(1))] );
loglog( dofDiagPlot(2:end,3,1), pathH1H2RelDistDiag(:,3,1),':xk','LineWidth',1,...
            'Displayname',['$n_\theta$ = ',num2str(nSList(3)),', $n_t$ = ',num2str(nTList(1))] );
loglog( dofDiagPlot(2:end,1,2), pathH1H2RelDistDiag(:,1,2),'-ok','LineWidth',1,...
            'Displayname',['$n_\theta$ = ',num2str(nSList(1)),', $n_t$ = ',num2str(nTList(2))] );
loglog( dofDiagPlot(2:end,2,2), pathH1H2RelDistDiag(:,2,2),'--ok','LineWidth',1,...
            'Displayname',['$n_\theta$ = ',num2str(nSList(2)),', $n_t$ = ',num2str(nTList(2))] );
loglog( dofDiagPlot(2:end,3,2), pathH1H2RelDistDiag(:,3,2),':ok','LineWidth',1,...
            'Displayname',['$n_\theta$ = ',num2str(nSList(3)),', $n_t$ = ',num2str(nTList(2))] );
hold off


% title('Relative H1H2 distance');
% ylabel('$||d_i - d_{i-1}||_{H1H2}/||d_i||_{H1H2}$', 'Interpreter','latex')
% xlabel('DOF = N*Nt')
hleg = legend('show','location','northeast')
set(hleg,'Interpreter','latex','FontSize',12)
set(gcf, 'Color','white');
% export_fig ConvergenceCircle2Wrap/pdf/RelH1H2DistDiagonal.pdf  
export_fig([plotDir,'/RelL2L2DistDiagonalPropellers.pdf'])


