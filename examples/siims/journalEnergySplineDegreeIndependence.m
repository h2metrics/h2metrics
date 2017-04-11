%% Energy convergence experiment (fixed parameters for theta) 

nTList = 1:1:4;
NtList = 10:10:110;


subdir = [prefixDir, 'EnergySplineDegree']; %Subdirectory for results
if ~exist(subdir, 'dir')
    mkdir(subdir);
end
%name = 'QuadPoints'; %name for the result files




f0 = @(t) 2*pi/(9.4248)/1/100*[17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - 2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + 13*sin(3*t) - sin(4*t) - 5*sin(5*t),...
 -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + 13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + 19*sin(3*t) - 8*sin(4*t)];
f1 = @(t) 2*pi/(9.4248)/100*[150*cos(-t - 3*pi/4), 150*sin(-t - 3*pi/4)] ;


splineData = constructEmptySplineData;
splineData.dSpace = 2;
splineData.quadDegree = [8, 6]; %Quadrature precission
splineData.N = 60; %Number of time control points
splineData.nS = 4; %time degree
splineData.a = [1, 0, 10^(-3)]; 

% Calculate Energies
for ii = 1:length(nTList)
    for jj =1:length(NtList)
                splineData.nT = nTList(ii);
                splineData.Nt = NtList(jj);
                splineData = constructKnots(splineData);
                splineData = setupQuadData(splineData);
                d0 = constructSplineApproximation(f0,splineData);
                d1 = constructSplineApproximation(f1,splineData);
                dInit = nonlinearPath(d0,d1,splineData);
                disp(['Case: nT = ',num2str(nTList(ii)),' ,Nt = ', num2str(NtList(jj))]);
                E = pathRiemH2Energy(dInit, splineData);
                filename =  strcat('nT',num2str(nTList(ii)),'Nt',num2str(NtList(jj)));
                savePath(filename,dInit,splineData,...
                    'workdir',subdir,'energy',E,'manualName',true);              
    end
end

%
energyArray = zeros(length(nTList),length(NtList),1);
for ii = 1:length(nTList)
    for jj =1:length(NtList)
        filename =   strcat('nT',num2str(nTList(ii)),'Nt',num2str(NtList(jj)));
        loaddir = [subdir,'/',filename,'.mat'];
        [dPathLoad, splineDataLoad,energyLoad] = loadPath(loaddir);
        %Save information
        energyArray(ii,jj,1) = energyLoad;
    end
end

diff = zeros(length(nTList),length(NtList)-1);
for ii=1:length(nTList)
    for jj=1:length(NtList)-1
    diff(ii,jj)= abs(energyArray(ii,jj+1,1)-energyArray(ii,jj,1))./energyArray(ii,jj,1);   
    end  
end

% Plot

% Plot parameters (in points)
lineWidth = 400;
figRelSize = 0.49;
figRatio = 4/3;
sx = figRelSize * lineWidth;
sy = sx / figRatio;

handle = figure( 'PaperUnits', 'points', 'PaperSize', [sx, sy], ...
                 'Units', 'points', 'Position', [0, 0, sx, sy], ...
                 'Color', 'white' );
handle.Visible = 'on';
%semilogy(diff(1,:),'k-o');
loglog(NtList(1:end-1),diff(1,:),'k-o');
hold on
loglog(NtList(1:end-1),diff(2,:),'k--s');
loglog(NtList(1:end-1),diff(3,:),'k:x');
loglog(NtList(1:end-1),diff(4,:),'k-.*');
set(gca, 'XDir')
hold off
set(gca, 'FontSize', 6);
AX=legend({'$n_{t}=1$', '$n_{t}=2$','$n_{t}=3$','$n_{t}=4$'},'Interpreter','latex','Location','southwest');
set(AX,'FontSize',7); 
export_fig([plotDir,'EnergyControllPointsFixedS.eps']);
close;






%% Nt FIXED
nSList = 2:1:5;
NsList = 10:10:110;

% Circle to wrap
f0 = @(t) 2*pi/(9.4248)/1/100*[17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - 2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + 13*sin(3*t) - sin(4*t) - 5*sin(5*t),...
 -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + 13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + 19*sin(3*t) - 8*sin(4*t)];
f1 = @(t) 2*pi/(9.4248)/100*[150*cos(-t - 3*pi/4), 150*sin(-t - 3*pi/4)] ;

splineData = constructEmptySplineData;
splineData.dSpace = 2;
splineData.quadDegree = [8, 6]; %Quadrature precission
splineData.Nt = 20 + 2; %Number of time control points
splineData.nT = 3; %time degree
splineData.a=[1, 0, 10^(-3)];

% Calculate Energies
for ii = 1:length(nSList)
    for jj =1:length(NsList)
                splineData.nS = nSList(ii);
                splineData.N = NsList(jj);
                splineData = constructKnots(splineData);
                splineData = setupQuadData(splineData);
                d0 = constructSplineApproximation(f0,splineData);
                d1 = constructSplineApproximation(f1,splineData);
                dInit = nonlinearPath(d0,d1,splineData);
                disp(['Case: nS = ',num2str(nSList(ii)),' ,Ns = ', num2str(NsList(jj))]);
                E = pathRiemH2Energy(dInit, splineData);
                filename =  strcat('nS',num2str(nSList(ii)),'Ns',num2str(NsList(jj)));
                savePath(filename,dInit,splineData,...
                    'workdir',subdir,'energy',E,'manualName',true);              
    end
end

energyArray = zeros(length(nSList),length(NsList),1);
for ii = 1:length(nSList)
    for jj =1:length(NsList)
        filename =   strcat('nS',num2str(nSList(ii)),'Ns',num2str(NsList(jj)));
        loaddir = [subdir,'/',filename,'.mat'];
        [dPathLoad, splineDataLoad,energyLoad] = loadPath(loaddir);
        %Save information
        energyArray(ii,jj,1) = energyLoad;
    end
end


diff = zeros(length(nSList),length(NsList)-1);
for ii=1:length(nSList)
    for jj=1:length(NsList)-1
    diff(ii,jj)= abs(energyArray(ii,jj+1,1)-energyArray(ii,jj,1))./energyArray(ii,jj,1);   
    end  
end

% Plot

% Plot parameters (in points)
lineWidth = 400;
figRelSize = 0.49;
figRatio = 4/3;
sx = figRelSize * lineWidth;
sy = sx / figRatio;

handle = figure( 'PaperUnits', 'points', 'PaperSize', [sx, sy], ...
                 'Units', 'points', 'Position', [0, 0, sx, sy], ...
                 'Color', 'white' );
handle.Visible = 'on';
%semilogy(diff(1,:),'k-o');
loglog(NsList(1:end-1),diff(1,:),'k-o');
hold on
loglog(NsList(1:end-1),diff(2,:),'k--s');
loglog(NsList(1:end-1),diff(3,:),'k:x');
loglog(NsList(1:end-1),diff(4,:),'k-.*');
set(gca, 'XDir')
hold off
AX=legend({'$n_{\theta}=2$', '$n_{\theta}=3$','$n_{\theta}=4$','$n_{\theta}=5$'},'Interpreter','latex','Location','southwest');
set(AX,'FontSize',7); 
set(gca, 'FontSize', 6);
export_fig([plotDir,'EnergyControllPointsFixedT.eps']);
close;





























    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
