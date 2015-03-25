%% Initialize
clear all;
close all;
for j=1:1


%% Define parameters
stopBeforeAMPL = false;
saveresults = true;
task = 'circle2wrap';

N = 10; %number of spacial control points, must be bigger than n+1
Nt_inner = 50; % Number of time control points not including the first and last
Nt = Nt_inner + 2; %Number of time control points
nS = 4; %spacial degree of spline
nT = 4; %time degree of spline
quadDegree = [8,4]; %Quadrature precission
s = 300; %Interpolationpoints

runfile='H2.run';
modfile='H2.mod';
datfile='H2.dat';
tabfile='H2.tab'; % no special characters allowed in this file name
tabledef={...
  'd[t,p,1]', 'dx';
  'd[t,p,2]', 'dy';
};
amploptions={...
  'option ipopt_options "max_iter=1000"'
};



%% Create or load the initial curves c0, c1

disp(['main.m, task =' task]);

switch task
    
case 'circle2circle' 
interpolS = linspace( 0, 2*pi, s+1); 
interpolS = interpolS(1:end-1);    
f0 = @(t) [cos(t);sin(t)]
f1 = @(t) 2*[cos(t);sin(t)] ;
c0=f0(interpolS); 
c1= f1(interpolS);


case 'circle2transcircle' 
interpolS = linspace( 0, 2*pi, s+1); 
interpolS = interpolS(1:end-1);     
f0 = @(t) [cos(t);sin(t)];
f1 = @(t) [cos(t)+5;sin(t)+5] ;
c0=f0(interpolS); 
c1= f1(interpolS);


case 'circle2wrap' 
interpolS = linspace( 0, 2*pi, s+1); 
interpolS = interpolS(1:end-1);     
wrap = @(t) 1/100*[17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - 2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + 13*sin(3*t) - sin(4*t) - 5*sin(5*t);...
 -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + 13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + 19*sin(3*t) - 8*sin(4*t)];
wrapcircle = @(t) 1/100*[150*cos(-t - 3*pi/4); 150*sin(-t - 3*pi/4)] ;
c0=wrap(interpolS); c1= wrapcircle(interpolS);


case 'cow2cat' 
interpolS = linspace( 0, 2*pi, 300+1); 
interpolS = interpolS(1:end-1);  
load('../02_Data/cow.txt');
load('../02_Data/cat.txt');
c0 = cow;
c1 = cat;

end



%% Set up spline and quadrature data 

disp('main.m, calling setupQuadData.m');
tic
splineData = struct( 'nS',nS,'nT',nT, 'N', N, 'Nt_inner',Nt_inner,...
    'quadDegree',quadDegree,'knotsS',[],'knotsT',[],'innerKnotsS',[],'innerKnotsT',[]);
splineData.knotsS =  [(-nS):(N+nS)]/N*2*pi; %normalize, domain of definition is [0,2*pi]x[0,2*pi]
%splineData.knotsT = [ zeros(nT,1), linspace(0,1,Nt_inner+2 - nT + 1), ones(nT,1)];
splineData.knotsT = [ zeros(1,nT), linspace(0,1,Nt- nT + 1), ones(1,nT)];
splineData.innerKnotsS = splineData.knotsS(nS+1:end-nS);
splineData.innerKnotsT = splineData.knotsT(nT+1:end-nT);
splineData.Nt_inner = Nt_inner;
[innerKnotsS, innerKnotsT, quadData] = setupQuadData(splineData);
toc

% Calculate the control points of the initial and final curves
B_interpol = spcol( splineData.knotsS, nS+1, brk2knt( interpolS, 1 ),'sparse');
B_interpol_p = [B_interpol(:,1:nS) + B_interpol(:,end-nS+1:end), B_interpol(:,nS+1:end-nS)];
d0 = B_interpol_p\c0';
d1 = B_interpol_p\c1';  

%% Write .dat file 

disp(['main.m, calling writedatfile.m, datfile = ' datfile]);
tic
writedatfile(d0,d1,quadData,datfile);
toc

%% Write .run file

disp(['main.m, calling writerunfile.m, runfile = ' runfile]);
tic
writerunfile(runfile, modfile, datfile, tabfile, tabledef, amploptions);
toc

%% Call AMPL
if stopBeforeAMPL 
    disp('main.m, stopped before running AMPL');
    break
end

disp('main.m, calling AMPL');
setenv('LD_LIBRARY_PATH','/users/bauerm/work/15_metrics_horeqnor/numerics/Ipopt/build/lib')
setenv('LD_RUN_PATH', '/users/bauerm/work/15_metrics_horeqnor/numerics/Ipopt/build/lib')
setenv('PATH', [getenv('PATH') ':/users/herman/COCONUT/bin:/users/bauerm/work/15_metrics_horeqnor/numerics/Ipopt/build/bin'])

% export LD_LIBRARY_PATH=/users/bauerm/work/15_metrics_horeqnor/numerics/Ipopt/build/lib
% export LD_RUN_PATH=/users/bauerm/work/15_metrics_horeqnor/numerics/Ipopt/build/lib
% export PATH=$PATH:/users/herman/COCONUT/bin:/users/bauerm/work/15_metrics_horeqnor/numerics/Ipopt/build/bin

tic  
[status, cmdout] = system(['ampl ' runfile],'-echo');
toc

if status ~= 0 
  error('error executing AMPL');
  return;
end

%% Read .tab file, extract spline data

disp(['main.m, calling writetabfile.m, tabfile = ' tabfile]);
tic
[data, variables] = readtabfile(tabfile);
toc
d=[squeeze(data(1,:))' squeeze(data(2,:))']; % size(d) = [nKP*nKT,2]




%% Plot snapshots of the curve
noPlotPointsS = 100;
plotPointsS = linspace(0,2*pi,noPlotPointsS);
noSnapshots = 9;
plotSnapShots = linspace(0,1,noSnapshots);
plotDataSH = setupPlotData( plotPointsS,plotSnapShots, splineData);
C_SH = plotDataSH.B*d;
figure(3)
for ii = 1:noSnapshots
    subplot( 3,3, ii)
    plot( C_SH( 1 + (ii-1)*plotDataSH.noPlotPointsS: (ii-1)*plotDataSH.noPlotPointsS + plotDataSH.noPlotPointsS, 1),...
        C_SH( 1 + (ii-1)*plotDataSH.noPlotPointsS: (ii-1)*plotDataSH.noPlotPointsS + plotDataSH.noPlotPointsS, 2) );
    axis equal
end


figure(2)
C_optimal = plotDataSH.B*d;
clf
hold on
for ii = 1:plotDataSH.noPlotPointsT
    plot( C_optimal( 1 + (ii-1)*plotDataSH.noPlotPointsS: (ii-1)*plotDataSH.noPlotPointsS + plotDataSH.noPlotPointsS, 1),...
        C_optimal( 1 + (ii-1)*plotDataSH.noPlotPointsS: (ii-1)*plotDataSH.noPlotPointsS + plotDataSH.noPlotPointsS, 2) );
    axis equal
end
hold off



%% Plot the speed of the path

speed = evaluateSpeed(d, quadData);
figure(4)
plot( quadData.points{2}', speed')
legend('Energy of the optimal path')

%% Save the results in a struct
%Save the results in a struct (overwrites on name conflict)
%Create a dynamic name for the results
if saveresults
    resultname = [task,'_','n',num2str(nS),'_','nT',num2str(nT),'_','N',num2str(N),'_','Nt',num2str(Nt),...
       'q',num2str(quadDegree(1)),num2str(quadDegree(2))];
    %Create an optimization result structure
    resultstruct = struct('d',d,'N',N,'Nt',Nt,'nS',nS,'nT',nT);
    %Save these in structure containing all results
    %Check if a file already exist, if so, append results
    if exist('spline_results.mat','file')
        load('spline_results.mat','spline_results')
        spline_results.(resultname) = resultstruct;
        save('spline_results.mat','spline_results');
    else
        spline_results.(resultname) = resultstruct;
        save('spline_results.mat','spline_results');
    end
end
clear all;
close all;

end
















