%% geodesicBvp
%
% Calculates the minimal geodesic betweeen 
%     d0 and d1 o G
% where d0, d1 are curves and G is a group of transformations. At the
% moment we support G to be any combination of
%   - Translations
%   - Rotations
%   - Reparametrizations (full diffeomorphism group)
%   - Rigid shifts of the parametrization (S^1 as a subgroup of Diff(S^1))
% 
% Input
%   d0, d1
%       Initial and final curves. Matrix of dimensions [N, dSpace].
%   splineData
%       General information about the splines used.
%   quadData, quadDataTensor
%       Precomputed spline collocation matrices at quadrature points.
%
% Optional parameters
%   options
%       Struct containing optimization options. Will be passed on to
%       subroutines. Contains the following fields:
%           useAmpl = {true, false (default)}
%           optTra = {true, false (default)}
%           optRot = {true, false (default)}
%           optDiff = {true, false (default)}
%           optShift = {true, false (default)}
%           maxIter = integer ([] for default value)
%           display = string
%               'off' for no output
%               '' for default of optimization routine
%               Any other string will be passed on
%   initPath
%       Guess for initial path.
%   gaInit
%       Guess for initial transformation of d1.
%
% Output
%   optE
%       Energy of the optimal path
%   optPath
%       Optimal path between d0 and d1 o optGa
%   optGa
%       Transformation between d0 and endpoint of optPath
%   info
%       Structure containing information about the minimization process
%
function [optE, optPath, optGa, info] = geodesicBvp(d0, d1, ...
    splineData, quadData, quadDataTensor, varargin)

%TODO: Add all decision functionality
optTra = false;
optRot = false;
optDiff = false;
optShift = false;

% Some code for handling optional inputs
ii = 1;
while ii <= length(varargin)
    if strcmpi(varargin{ii}, 'options')
        ii = ii + 1;
        options = varargin{ii};
    end
    ii = ii + 1;  
end

% Enforce default options
if ~isfield(options, 'useAmpl')
    options.useAmpl = false;
end

if ~options.useAmpl
    %TODO: Create real geodesicBvpMatlab functionality
%     [optE, optPath, optGa, info] = geodesicBvpMatlab(d0, d1, ...
%         splineData, quadData, quadDataTensor, varargin);
    
    %Do basic parametrized energy
    if ((~optTra) && (~optRot) && (~optDiff) && (~optShift)) 
        [optE, optPath, optGa, info] = geodesicBvpFminunc(d0, d1, ...
            splineData, quadData, quadDataTensor, varargin);
    end
else
    error('Ampl version not implemented yet...');
end

end

%Actual optimization function
function [E_geo, dPath, opt_grad,opt_output] = geodesicBvpFminunc(d0,d1,...
        splineData,quadData,quadDataTensor,varargin);
% Compute the minimal geodesic connecting the splines given by d0 and d1.%
%
% Input: 
%       d0, [NxdSpace], first set of control points
%       d1, [NxdSpace], second set of control points
%       splineData,
%       quadData,
%       quadDataTensor
%       'option' = 'value', list of options on the form 
%       
% Output:
%       E_geo, optimal energy
%       dPath, optimal path
%
d_linear = linearPath(d0,d1,splineData);
d_init = d_linear( splineData.N+1:end-splineData.N,:);

%Some code for handling optional inputs
% ii = 1;
% while ii <= length(varargin)-1
%     if (isa(varargin{ii},'char') && isa(varargin{ii+1},'char'))
%         switch (lower(varargin{ii}))
%             case 'display' %Don't return gradient terms from end curves
%                 options = optimoptions(options,'Display',lower(varargin{ii+1}));
%             case 'plotfval' %value = 'true'/'false'
%                 if strcmpi(varargin{ii+1},'true')
%                     options = optimoptions(options,'PlotFcns', @optimplotfval);
%                 end
%             case 'algorithm'
%                 options = optimoptions(options,'Algorithm',lower(varargin{ii+1}));
%             case 'gradobj'
%                 options = optimoptions(options,'GradObj', lower(varargin{ii+1}));
%         end
%     elseif (isa(varargin{ii},'char') && isa(varargin{ii+1},'numeric'))
%         switch (lower(varargin{ii}))
%             case 'tolfun'
%                 options = optimoptions(options,'TolFun', varargin{ii+1});
%             case 'tolx'
%                 options = optimoptions(options,'TolX', varargin{ii+1});
%             case 'init'
%                 d_init = varargin{ii+1};
%             case 'maxfunevals'
%                 options = optimoptions(options, 'MaxFunEvals',varargin{ii+1});
%         end
%     end
%     ii = ii + 1;
% end

options = optimoptions('fminunc');
options = optimoptions(options,'GradObj', 'on');
options = optimoptions(options,'Hessian', 'on');
options = optimoptions(options,'HessMult',[]);
options = optimoptions(options,'Algorithm', 'trust-region');
options = optimoptions(options,'Display','iter');
% options = optimoptions(options,'MaxPCGIter',4000);
% options = optimoptions(options,'TolFun',1e-6);


[d_optimal,E_geo,exitflag,output,grad] = ...
fminunc(@(d_opt)energyH2([d0;d_opt;d1],splineData,quadDataTensor),...
    d_init,options);

dPath = [d0;d_optimal;d1];

if nargout > 2 %provide output options
    opt_grad = grad;
end
if nargout > 3 %provide output options
    opt_output = output;
end

end

