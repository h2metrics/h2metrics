function [E_geo, dPath] = geodesicBVP_ampl(d0,d1,splineData,quadData,quadDataTensor);
% Compute the minimal geodesic connecting the splines given by d0 and d1 using AMPL to solve the minimization problem.%
%
% Input: 
%       d0, [NxdSpace], first set of control points
%       d1, [NxdSpace], second set of control points
%       splineData,
%       quadData,
%       quadDataTensor
%
%       
% Output:
%       E_geo, optimal energy
%       dPath, optimal path
%
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

%% Write .dat file
disp(['main.m, calling writedatfile.m, datfile = ' datfile]);
tic
writedatfile(d0,d1,splineData,quadData,quadDataTensor,datfile);
toc


%% Write .run file

disp(['geodesicBVP_ampl.m, calling writerunfile.m, runfile = ' runfile]);
tic
writerunfile(runfile, modfile, datfile, tabfile, tabledef, amploptions);
toc

disp('geodesicBVP_ampl.m, calling AMPL');
setenv('LD_LIBRARY_PATH','/users/bauerm/work/15_metrics_horeqnor/numerics/Ipopt/build/lib')
setenv('LD_RUN_PATH', '/users/bauerm/work/15_metrics_horeqnor/numerics/Ipopt/build/lib')
setenv('PATH', [getenv('PATH') ':/users/herman/COCONUT/bin:/users/bauerm/work/15_metrics_horeqnor/numerics/Ipopt/build/bin'])

tic  
[status, cmdout] = system(['ampl ' runfile],'-echo');
toc

if status ~= 0 
  error('error executing AMPL');
  return;
end

%% Read .tab file, extract spline data

disp(['geodesicBVP_ampl.m, calling writetabfile.m, tabfile = ' tabfile]);
tic
[data, variables] = readtabfile(tabfile);
toc
dPath=[squeeze(data(1,:))' squeeze(data(2,:))']; 
E_geo = energyH2(dPath,splineData,quadDataTensor)


end

