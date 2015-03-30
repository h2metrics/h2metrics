%% Initialize
clear all;
close all;


%% Define parameters

stopBeforeAMPL = false;
saveresults = true;

task = struct(...
    'shapes','circle2wrap',...          % choice of initial and final curve
    'pHor', 0,...                       % pHor=1 means horizontal energy
    'pH0', 1,...                        % the coefficient of the H^0 term of the metric
    'pKappa', 0,...                     % the coefficient of the kappa^2 term of the metric
    'pH2',1,...                         % the coefficient of the H^2 term of the metric
    'nKU',16,...                        % number of spacial control points, must be bigger than n+1
    'nKT',11,...                        % number of time control points 
    'dU',4,...                          % spacial degree of spline
    'dT',4,...                          % time degree of spline
    'quadDegree',[8,4],...              % quadrature degrees in space and time
    'nrPlotPoints',[100,9],...          % number of points in space and time used for plotting results
    'nrInterpolationPoints',300, ...    % used to generate initial and final shapes
    'runfile','H2.run',...              % file name of the .run file, including the file extension
    'modfile','H2.mod',...              % file name of the .mod file, including the file extension
    'datfile','H2.dat',...              % file name of the .dat file, including the file extension
    'tabfile','H2.tab',...              % file name of the .tab file, including the file extension
    'resultsfile','spline_results.mat',... % file name of the .mat file used to store matlab results
    'tabledef',[],...                   % definition of the table in the .tab file
    'amploptions',[],...                % options to be written in the .run file
    'splineData',[],...                 % collocation matrices and quadrature weights
    'splineDataSH',[],...               % the same, but for plotting 'snapshots' of the path
    'results',[]...
);

task.tabledef = { ...
    'd[t,u,1]','dx'; ...
    'd[t,u,2]','dy'
};

task.amploptions = { ...
    'option ipopt_options "max_iter=10"'
};

task.results = struct(...
    'c0',   [],...                      % initial curve
    'c1',   [],...                      % final curve
    'd0',   [],...                      % initial controls
    'd1',   [],...                      % final controls
    'd',   [] ...                       % all controls (the optimal ones)
);

%% Set up spline and quadrature data 

disp('main.m, setting up spline and quadrature data');
task.splineData = setupSplineData(task);

%% Create or load the initial curves c0, c1

disp(['main.m, loading shapes ' task.shapes]);

switch task.shapes
    
case 'circle2circle' 
x = linspace( 0, 2*pi, task.nrInterpolationPoints+1); 
x = x(1:end-1);    
f0 = @(t) [cos(t);sin(t)]
f1 = @(t) 2*[cos(t);sin(t)] ;
task.results.c0=f0(x); 
task.results.c1= f1(x);
clear f0 f1;


case 'circle2transcircle' 
x = linspace( 0, 2*pi, task.nrInterpolationPoints+1); 
x = x(1:end-1);     
f0 = @(t) [cos(t);sin(t)];
f1 = @(t) [cos(t)+5;sin(t)+5] ;
task.results.c0=f0(x); 
task.results.c1= f1(x);
clear f0 f1;

case 'circle2wrap' 
x = linspace( 0, 2*pi, task.nrInterpolationPoints+1); 
x = x(1:end-1);     
wrap = @(t) 1/100*[17 - 34*cos(t) - 79*cos(2*t) + 7*cos(3*t) - 2*cos(4*t) + 3*cos(5*t) - 24*sin(t) - 86*sin(2*t) + 13*sin(3*t) - sin(4*t) - 5*sin(5*t);...
 -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + 13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + 19*sin(3*t) - 8*sin(4*t)];
wrapcircle = @(t) 1/100*[150*cos(-t - 3*pi/4); 150*sin(-t - 3*pi/4)] ;
task.results.c0=wrap(x); 
task.results.c1= wrapcircle(x);
clear wrap wrapcircle


case 'cow2cat' 
x = linspace( 0, 2*pi, 300+1); 
x = x(1:end-1);  
load('../02_Data/cow.txt');
load('../02_Data/cat.txt');
task.results.c0 = cow;
task.results.c1 = cat;

end

% Calculate the control points of the initial and final curves
B_interpol = spcol( task.splineData.KU, task.dU+1, brk2knt( x, 1 ),'sparse');
B_interpol_p = [B_interpol(:,1:task.dU) + B_interpol(:,end-task.dU+1:end), B_interpol(:,task.dU+1:end-task.dU)];
task.results.d0 = B_interpol_p\task.results.c0';
task.results.d1 = B_interpol_p\task.results.c1';  

clear x B_interpol B_interpol_p 

%% Write .dat file 

disp(['main.m, calling writedatfile.m, datfile = ' task.datfile]);
tic
writedatfile(task);
toc

%% Write .run file

disp(['main.m, calling writerunfile.m, runfile = ' task.runfile]);
tic
writerunfile(task);
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
[status, cmdout] = system(['ampl ' task.runfile],'-echo');
toc

if status ~= 0 
  error('error executing AMPL');
  return;
end

%% Read .tab file, extract spline data

disp(['main.m, calling writetabfile.m, tabfile =' task.tabfile]);
tic
[data, variables] = readtabfile(task.tabfile);
toc
task.results.d=[squeeze(data(1,:))' squeeze(data(2,:))']; % size(task.d) = [nKP*nKT,2]




%% Plot snapshots of the curve

% call setupSplineDataSH;
task.splineDataSH = setupSplineDataSH(task);

C_SH = task.splineDataSH.B*task.results.d;
nQU = task.nrPlotPoints(1); nQT = task.nrPlotPoints(2);
figure(3)
for ii = 1:nQT
    subplot( 3,3, ii)
    plot( C_SH( 1 + (ii-1)*nQU: (ii-1)*nQU + nQU, 1),...
        C_SH( 1 + (ii-1)*nQU: (ii-1)*nQU + nQU, 2) );
    axis equal
end

figure(2)
clf
hold on
for ii = 1:task.splineDataSH.nQT
    plot( C_SH( 1 + (ii-1)*nQU: (ii-1)*nQU + nQU, 1),...
        C_SH( 1 + (ii-1)*nQU: (ii-1)*nQU + nQU, 2) );
    axis equal
end
hold off



%% Plot the speed of the path

speed = evaluateSpeed(task);
figure(4)
plot( task.splineData.QT', speed')
legend('Energy of the optimal path')

%% Save the results in a struct
%Save the results in a struct (overwrites on name conflict)
%Create a dynamic name for the results
if saveresults
    resultname = [task.shapes,'_','n',num2str(task.dU),'_','dT',num2str(task.dT),'_','N',num2str(task.nKU),'_','Nt',num2str(task.nKT),...
       'q',num2str(task.quadDegree(1)),num2str(task.quadDegree(2))];
    %Save these in structure containing all results
    %Check if a file already exist, if so, append results
    if exist(task.resultsfile,'file')
        load(task.resultsfile,'spline_results')
        spline_results.(resultname) = task;
        save(task.resultsfile,'spline_results');
    else
        spline_results.(resultname) = task;
        save(task.resultsfile,'spline_results');
    end
end
clear all;
close all;

















