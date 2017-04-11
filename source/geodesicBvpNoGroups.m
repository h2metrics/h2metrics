%% geodesicBvpNoGroups
%
% Calculates the minimal geodesic betweeen 
%     d0 and d1 o G
% where d0, d1 are curves and G is a group of transformations. This
% function supports G to be any combination of
%   - Translations
%   - Rotations
%   - Rigid shifts of the parametrization
% 
% Input
%   d0, d1
%       Initial and final curves. Matrix of dimensions [N, dSpace].
%   splineData
%       General information about the splines used.
%
% Optional parameters
%   options
%       Struct containing optimization options. Uses the following fields:
%           optTra = false
%           optRot = false
%           optDiff = false
%           optShift = false
%           maxIter = integer ([] for default value)
%           display = string
%               'off' for no output
%               '' for default of optimization routine
%               Any other string will be passed on
%           usePrecond = {true, false (default)}
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
function [optE, optPath, optGa, info] = geodesicBvpNoGroups(d0, d1, ...
    splineData, varargin)

%% Default parameters
usePrecond = false;

options = [];
dInitPath = [];

%% Extract parameters
N = splineData.N;
Nt = splineData.Nt;
dSpace = splineData.dSpace;

%% Some code for handling optional inputs
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'options'
                ii = ii + 1;
                options = varargin{ii};
            case 'initpath'
                ii = ii + 1;
                dInitPath = varargin{ii};
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1;  
    end
end

%% Set options
if isfield(options, 'usePrecond')
    usePrecond = options.usePrecond;
end
   
minOptions = optimoptions('fminunc');
minOptions = optimoptions(minOptions,'Algorithm', 'trust-region');
minOptions = optimoptions(minOptions,'DerivativeCheck', 'off');
% minOptions = optimoptions(minOptions,'PlotFcns', @optimplotfval);
minOptions = optimoptions(minOptions,'GradObj', 'on');
minOptions = optimoptions(minOptions,'Hessian', 'on');
minOptions = optimoptions(minOptions,'MaxFunEvals', 1000000);
if isfield(options, 'display')
    minOptions = optimoptions(minOptions, 'Display', options.display);
end
if isfield(options, 'tolFun')
    minOptions = optimoptions(minOptions,'TolFun', options.tolFun);
end
if isfield(options, 'tolX')
    minOptions = optimoptions(minOptions,'TolX', options.tolX);
end
if isfield(options, 'maxIter')
    minOptions = optimoptions(minOptions, 'maxIter', options.maxIter);
end

%% Create initial guess for path if not provided one
if isempty(dInitPath)
    dInitPath = linearPath(d0, d1, splineData);
end

%% Preconditioning
if usePrecond
    % Scale parameter for preconditioner
    lambda = ( curveLength(d0, splineData)/(2*pi) + ...
               curveLength(d1, splineData)/(2*pi) ) / 2;
           
    [~, ~, ~, PinvBlock] = createPreconditioner(lambda, splineData);
else
    % Identity matrix as preconditioner
    PinvBlock = speye(N*(Nt-2)*dSpace);
end

%% Setup optimization
coeffInit = reshape( dInitPath(N+1:end-N, :), ...
                     [N*(Nt-2)*dSpace, 1] );
coeffInit = PinvBlock \ coeffInit; % Preconditioner

params = struct('d0', d0, 'd1', d1, 'PinvBlock', PinvBlock);

Fopt = @(coeffs) energyH2( coeffs, params, splineData );

problem = struct( 'objective', Fopt, 'x0', coeffInit, ...
                  'options', minOptions, 'solver', 'fminunc' );
                  
%% Optimize
% tic
[coeffOptimal, optE, exitflag, output] = fminunc( problem );
% toc

%% Save results
% Create transformation struct
optGa = struct( 'phi', [], 'beta', [], 'v', [], 'alpha', []);

coeffOptimal = PinvBlock * coeffOptimal;

optPath = [ d0; ...
          reshape(coeffOptimal(1:N*(Nt-2)*dSpace), [N*(Nt-2), dSpace]); ...
          d1 ];
      
info = struct( 'exitFlag', exitflag, ...
               'noIter', output.iterations );

end
