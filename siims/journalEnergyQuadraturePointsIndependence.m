%% Plot parameters (in points)
lineWidth = 400;
figRelSize = 0.49;
rigidOptions = struct( 'optDiff', false, ...
                       'optTra', false, ...
                       'optRot', false, ...
                       'optShift', true );




%%
subdir = 'QuadraturePointsIndependence'; %Subdirectory for results
name = 'QuadPoints'; %name for the result files
quadDegreeSList = 3:1:18;
quadDegreeTList = 1:1:16;
%Number of curves
noCu=4;


% Circle to wrap
f0 = @(t) 2*pi/(9.4248)/1/100*[17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - 2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + 13*sin(3*t) - sin(4*t) - 5*sin(5*t),...
 -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + 13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + 19*sin(3*t) - 8*sin(4*t)];
f1 = @(t) 2*pi/(9.4248)/100*[150*cos(-t - 3*pi/4), 150*sin(-t - 3*pi/4)] ;

splineData = constructEmptySplineData;
splineData.dSpace = 2;
splineData.N = 60; %no. control points, must be bigger than n+1
splineData.nS = 4; %spacial degree
splineData.Nt = 20 + 2; %Number of time control points
splineData.nT = 3; %time degree



%% 
for ii = 1:length(quadDegreeSList)
                splineData.quadDegree = [quadDegreeSList(ii), quadDegreeTList(ii)]; %Quadrature precission
                splineData = constructKnots(splineData);
                splineData.a=[1, 0, 10^(-3)]; 
                [quadData, quadDataTensor,splineData] = setupQuadData(splineData);
                d0 = constructSplineApproximation(f0,splineData);
                d1 = constructSplineApproximation(f1,splineData);
                dInit = nonlinearPath(d0,d1,splineData);
                disp(['Case: quadDegreeSList = ',num2str(quadDegreeSList(ii)),' ,quadDegreeTList = ', num2str(quadDegreeTList(ii))]);
                E = pathRiemH2Energy(dInit, splineData,quadData, quadDataTensor);
                quadpoints = splineData.quadData.noQuadPointsS*splineData.quadData.noQuadPointsT;
                quadpoints= [splineData.quadData.noQuadPointsS,splineData.quadData.noQuadPointsT];
                filename =  strcat('Case',num2str(1),'quadDegreeS',num2str(quadDegreeSList(ii)),'quadDegreeT',num2str(quadDegreeTList(ii)));
                savePath2(filename,dInit,splineData,...
                    'workdir',subdir,'energy',E,'quadpoints',quadpoints,'manualName',true);              
end

for ii = 1:length(quadDegreeSList)
                splineData.quadDegree = [quadDegreeSList(ii), quadDegreeTList(ii)]; %Quadrature precission
                splineData = constructKnots(splineData);
                splineData.a=[1, 0, 10^(-3)]; 
                [quadData, quadDataTensor,splineData] = setupQuadData(splineData);
                dList = loadDataSet( 'basic', splineData,'fft' );
                d0 = dList{6};
                d1 = dList{7};
                dInit = nonlinearPath(d0,d1,splineData);
                disp(['Case: quadDegreeSList = ',num2str(quadDegreeSList(ii)),' ,quadDegreeTList = ', num2str(quadDegreeTList(ii))]);
                E = pathRiemH2Energy(dInit, splineData,quadData, quadDataTensor);
                quadpoints = splineData.quadData.noQuadPointsS*splineData.quadData.noQuadPointsT;
                filename =  strcat('Case',num2str(2),'quadDegreeS',num2str(quadDegreeSList(ii)),'quadDegreeT',num2str(quadDegreeTList(ii)));
                savePath2(filename,dInit,splineData,...
                    'workdir',subdir,'energy',E,'quadpoints',quadpoints,'manualName',true);              
end




for ii = 1:length(quadDegreeSList)
                splineData.quadDegree = [quadDegreeSList(ii), quadDegreeTList(ii)]; %Quadrature precission
                splineData = constructKnots(splineData);
                splineData.a=[1, 0, 10^(-5)]; 
                [quadData, quadDataTensor,splineData] = setupQuadData(splineData);
                dataDir = '../data/';
                dList = loadDataSet( 'corpus_callosum_tilak', splineData, dataDir, ...
                     'class', 'normal', 'ind', [1, 2], 'constSpeed');
                d0 = dList{1};
                d1 = dList{2};
                meanLength = ( curveLength(d1, splineData, quadData) + ...
                curveLength(d1, splineData, quadData) ) / 2;
                d0 = d0 / meanLength * 2*pi;
                d1 = d1 / meanLength * 2*pi;
                [dList, ~] = rigidAlignment( {d0, d1}, splineData, quadData, ...
                          'options', rigidOptions );
                d1 = dList{2};             
                dInit = nonlinearPath(d0,d1,splineData);
                disp(['Case: quadDegreeSList = ',num2str(quadDegreeSList(ii)),' ,quadDegreeTList = ', num2str(quadDegreeTList(ii))]);
                E = pathRiemH2Energy(dInit, splineData,quadData, quadDataTensor);
                quadpoints = splineData.quadData.noQuadPointsS*splineData.quadData.noQuadPointsT;
                filename =  strcat('Case',num2str(3),'quadDegreeS',num2str(quadDegreeSList(ii)),'quadDegreeT',num2str(quadDegreeTList(ii)));
                savePath2(filename,dInit,splineData,...
                    'workdir',subdir,'energy',E,'quadpoints',quadpoints,'manualName',true);              
end

%% Load results
energyArray = zeros(length(quadDegreeSList),3,noCu);
diff = zeros(length(quadDegreeSList)-1,noCu);
for jj=1:noCu
    for ii = 1:length(quadDegreeSList)
        filename =  strcat('Case',num2str(jj),'quadDegreeS',num2str(quadDegreeSList(ii)),'quadDegreeT',num2str(quadDegreeTList(ii)));
        loaddir = [subdir,'/',filename,'.mat'];
        [dPathLoad, splineDataLoad,energyLoad,quadpointsLoad] = loadPath2(loaddir);
        %Save information
        energyArray(ii,1,jj) = energyLoad;
        energyArray(ii,2,jj) = splineDataLoad.quadDegree(1);
        energyArray(ii,3,jj) = splineDataLoad.quadDegree(2);
        diff(:,jj)= abs(energyArray(2:end,1,jj)-energyArray(1:end-1,1,jj))./energyArray(1:end-1,1,jj);
    end
end    
diff2 = zeros((length(quadDegreeSList))/2-1,noCu);

for jj=1:(length(quadDegreeSList))/2-1
    diff2(jj,:)= diff(2*jj-1,:);
end

%% Plot
figRatio = 4/3;
sx = figRelSize * lineWidth;
sy = sx / figRatio;

handle = figure( 'PaperUnits', 'points', 'PaperSize', [sx, sy], ...
                 'Units', 'points', 'Position', [0, 0, sx, sy], ...
                 'Color', 'white' );
handle.Visible = 'on';
semilogy(diff2(:,1),'k-o');
hold on
semilogy(diff2(:,2),'k--s');
semilogy(diff2(:,3),'k:x');
%semilogy(diff2(:,4),'k-.*');
set(gca, 'XDir')
set(gca,'XTick',[1,3,5,7]);
set(gca, 'XTickLabel', {'(1,2)';'(3,4)';'(5,6)';'(7,8)'});
hold off
AX=legend({'Circle to wrap', '3-blade to 4-blade','Corpus callosum 1 to 2'},'Location','southwest');
set(AX,'FontSize',5); 
export_fig('quadraturePointIndependence.eps');
close;










