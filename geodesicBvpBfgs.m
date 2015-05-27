function [E_geo,dPath,varargout] = geodesicBvpBfgs(d0,d1,splineData,quadData,quadDataTensor,varargin)
% Compute the minimal geodesic connecting the splines given by d0 and d1 using AMPL to solve the minimization problem.%
%
% Input: 
%       d0, [NxdSpace], first set of control points
%       d1, [NxdSpace], second set of control points
%       splineData,
%       quadData,
%       quadDataTensor
%       varargin .......... Boole variable, that controlls wether datfile2
%       is written (value = 1 --> datfile exists). 
%       
%
%       
% Output:
%       E_geo, optimal energy
%       dPath, optimal path
%       varargout, exit code from AMPL

runfile='H2Bfgs.run';
datfile1='H2.dat';
datfile2='H2_tensor.dat';
tabfile='H2.tab'; % no special characters allowed in this file name
tabledef={...
      'd[t,p,1]', 'dx';
      'd[t,p,2]', 'dy';
              };
amploptions={...
      'option lbfgsb_options "pgtol=1e-04 maxit=1000";'};  
datfileexists = 0;
mintrans = 0;
minscale = 0;
minrot = 0;
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'datfileexists'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    datfileexists = logical(varargin{ii});
                else
                    error('Invalid value for option ''datfileexists''.');
                end
            case 'mintrans'
                ii = ii +1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    mintrans = logical(varargin{ii});
                else
                    error('Invalid value for option ''mintrans''.');
                end
            case 'minscale'
                ii = ii +1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    minscale = logical(varargin{ii});
                else
                    error('Invalid value for option ''minscale''.');
                end
            case 'minrot'
                ii = ii +1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    minrot = logical(varargin{ii});
                else
                    error('Invalid value for option ''minrot''.');
                end
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});    
        end
    end
    ii = ii + 1; 
end




if minscale == 1;
    modfile = 'H2Scale.mod';
    disp(['modfile H2scale.mod']);
else 
    modfile = 'H2.mod';
end    
            
persistent lasttimestamp;
if lasttimestamp == quadData.timestamp
    datfileexists = 1;
else
    lasttimestamp = quadData.timestamp;
    datfileexists = 0;
end
    

%% Write datfile1
disp(['geodesicBvpBfgs.m, calling writeDatFile1.m, datfile = ' datfile1]);
writeDatFile1(d0,d1,splineData,quadData,mintrans,minrot,datfile1);


%% Write datfile2
if datfileexists == 1;
    disp(['datfile2 exists']);
else 
    disp(['geodesicBvpBfgs.m, calling writeDatFile2.m, datfile = ' datfile2]);
    writeDatFile2(splineData,quadData,quadDataTensor,datfile2);
end
%% Write .run file
disp(['geodesicBvpBfgs.m, calling writeRunFile.m, runfile = ' runfile]);
writeRunFileBfgs(runfile, modfile, datfile1,datfile2, tabfile, tabledef, amploptions);

disp('geodesicBvpBfgs.m, calling AMPL');
%setenv('LD_LIBRARY_PATH','/users/bauerm/work/15_metrics_horeqnor/numerics/Ipopt/build/lib')
%setenv('LD_RUN_PATH', '/users/bauerm/work/15_metrics_horeqnor/numerics/Ipopt/build/lib')
%setenv('PATH', [getenv('PATH') ':/users/herman/COCONUT/bin:/users/bauerm/work/15_metrics_horeqnor/numerics/Ipopt/build/bin'])

tic  
[status, cmdout] = system(['ampl ' runfile],'-echo');
toc

if status ~= 0
    warning('geodesicBvpBfgs: Error executing AMPL');
else
    if ~isempty(strfind(cmdout,'L-BFGS-B: more than maxit iterations'))
        warning('geodesicBvpAmpl: Maximum Number of Iterations Exceeded');
        status = 1;
    else if isempty(strfind(cmdout,'L-BFGS-B: CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTO'))
        status = 1;
        warning('geodesicBvpAmpl: Unknown exit code from AMPL');
    end
end

if nargout>2 
    varargout{1} = status;
end

%% Read .tab file, extract spline data

disp(['geodesicBvpBfgs.m, calling readTabFile.m, tabfile = ' tabfile]);

[data, variables] = readTabFile(tabfile);

dPath=[squeeze(data(1,:))' squeeze(data(2,:))']; 
E_geo = pathRiemH2Energy( dPath, splineData,quadData, quadDataTensor);
%L_geo = pathRiemH2Length( dPath, splineData,quadData, quadDataTensor)
end

