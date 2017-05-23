%% geodesicBvp
%
% Minimizes the varifold--H^2 inexact boundary value problem
% between d0, d1\circ G are curves and G is a group of transformations using 
% the HANSO LM-BFGS method This
% function supports G to be any combination of
%   - Translations
%   - Rotations
%
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
%           optTra = {true, false (default)}
%           optRot = {true, false (default)}
%           maxIter = integer ([] for default value)
%           display = string
%               'off' for no output
%               '' for default of optimization routine
%               Any other string will be passed on
%   initPath
%       Guess for initial path. path needs to have the same splineData as splineData. 
%       Default initial path: constant path d0
%   gaInit
%       Guess for initial transformation of d1.
%   multigrid
%       Uses a lower resolution splineData for a preComputation
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
function [optE, optPath, optGa, info] = geodesicBvpVarifold(d0, d1, ...
    splineData, options, varargin)

%% Default parameters
dInitPath = [];
initGa = struct('beta', [], 'v', []);

% Default options
optRot = false;
optTra = false;

%% Read initial data
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'initpath'
                ii = ii + 1;
                dInitPath = varargin{ii};
            case 'initga'
                ii = ii + 1;
                initGa = varargin{ii};
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end  
    end
    ii = ii + 1;
end

% Set options
if isfield( options, 'optRot' )
    optRot = options.optRot;
end
if isfield( options, 'optTra' )
    optTra = options.optTra;
end

dSpace = splineData.dSpace;

%% Set options  
minOptions = optimoptions('fminunc');
minOptions = optimoptions(minOptions,'Algorithm', 'quasi-newton');
minOptions = optimoptions(minOptions,'DerivativeCheck', 'off');
minOptions = optimoptions(minOptions,'PlotFcns', @optimplotfval);
minOptions = optimoptions(minOptions,'GradObj', 'on');
minOptions = optimoptions(minOptions,'Hessian', 'off');
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

%% Create initial guess for path and gamma if not provided one
if isempty(dInitPath)
    dInitPath = linearPath(d0, d0, splineData);
end
if isfield(initGa, 'beta') && ~isempty(initGa.beta)
    initBeta = initGa.beta;
else
    initBeta = 0;
end
if isfield(initGa, 'v') && ~isempty(initGa.v)
    initV = initGa.v;
else
    initV = [0; 0];
end

coeffInit = [reshape(dInitPath(splineData.N+1:end,:),[],1); initBeta; initV];

%% Setup optimization
Fopt = @(coeff) energyH2Varifold( ...
    [d0; reshape(coeff(1:end-3), [],2)],d1,coeff(end-2),coeff(end-1:end),...
    splineData,'optRot',optRot,'optTra',optTra);

problem = struct( 'objective', Fopt, 'x0', coeffInit, ...
                  'options', minOptions, 'solver', 'fminunc' );
                            
%% Optimize
% tic
% [coeffOptimal, optE, exitflag, output] = fminunc( problem );
% toc

%% Setup HANSO

% FoptHANSO = @(coeff,pars) Fopt(coeff);

pars = struct(); optionsHANSO = struct();
pars.nvar = length(coeffInit);
pars.fgname = @energyH2VarifoldHANSO; %[f,g] = fgtest(x,pars)
pars.splineData = splineData;
pars.optRot = optRot;
pars.optTra = optTra;
pars.d0 = d0;
pars.d1 = d1;


optionsHANSO.x0 = coeffInit;
optionsHANSO.normtol = 1e-3;
optionsHANSO.maxit = 1000;
optionsHANSO.nvec = 500; %0 is full bfgs
optionsHANSO.fvalquit = 0;
optionsHANSO.prtlevel = 1; % also 0,2 

%% Optimize HANSO

[coeffOptimal, optE, infoHanso] = hanso(pars, optionsHANSO);

output.iterations = 1000; % ?????

%% Create output
% Transformation struct
optGa = struct( 'beta', [], 'v', [] );

if optTra
    optGa.v = coeffOptimal(end-dSpace+1:end);
end
if optRot
    optGa.beta = coeffOptimal(end-dSpace);
end

optPath = [ d0; reshape( coeffOptimal(1:end-dSpace-1), [], 2) ];
           
info = struct( 'infoHanso', infoHanso, ...
               'noIter', output.iterations ); 

end
