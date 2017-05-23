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
initPath = [];
initGa = struct('beta', [], 'v', []);

dSpace = splineData.dSpace;

%% Read initial data
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'initpath'
                ii = ii + 1;
                initPath = varargin{ii};
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
else
    optRot = false;
end
if isfield( options, 'optTra' )
    optTra = options.optTra;
else
    optTra = false;
end

%% Create initial guess for path if not provided one
if isempty(initPath)
    initPath = linearPath(d0, d0, splineData); % This is constant path
end
if isfield(initGa, 'beta') && ~isempty(initGa.beta)
    initBeta = initGa.beta;
else
    initBeta = 0;
end
if isfield(initGa, 'v') && ~isempty(initGa.v)
    initV = initGa.v;
else
    initV = zeros(dSpace, 1);
end

coeffInit = [ reshape(initPath(splineData.N+1:end, :), [], 1); ...
              initBeta; ...
              initV ];

%% Setup HANSO
Fopt = @(coeff, pars) energyH2Varifold( ...
    [pars.d0; reshape(coeff(1:end-3), [],2)], pars.d1, ...
    coeff(end-dSpace), coeff(end-dSpace+1:end), ...
    pars.splineData, 'optRot', pars.optRot, 'optTra', pars.optTra);

pars = struct();
pars.nvar = length(coeffInit);
pars.fgname = Fopt; %[f,df] = fgtest(x,pars)
pars.splineData = splineData;
pars.optRot = optRot;
pars.optTra = optTra;
pars.d0 = d0;
pars.d1 = d1;

optionsHANSO = struct();
optionsHANSO.x0 = coeffInit;
if isfield(options, 'hansoNormTol')
    optionsHANSO.normtol = options.hansoNormTol;
else
    optionsHANSO.normtol = 1e-3;
end
if isfield(options, 'hansoMaxIt')
    optionsHANSO.maxit = options.hansoMaxIt;
else
    optionsHANSO.maxit = 1000;
end
if isfield(options, 'hansoNvec')
    optionsHANSO.nvec = options.hansoNvec;
else
    optionsHANSO.nvec = 500; % 0 is full bfgs
end
optionsHANSO.fvalquit = 0;
if isfield(options, 'hansoPrtLevel')
    optionsHANSO.prtlevel = options.hansoPrtLevel;
else
    optionsHANSO.prtlevel = 1; % also 0, 2
end

%% Call HANSO
[optCoeff, optE, infoHanso] = hanso(pars, optionsHANSO);

%% Create output
% Transformation struct
optGa = struct( 'beta', [], 'v', [] );

if optTra
    optGa.v = optCoeff(end-dSpace+1:end);
end
if optRot
    optGa.beta = optCoeff(end-dSpace);
end

optPath = [ d0; reshape( optCoeff(1:end-dSpace-1), [], 2) ];
           
info = struct( 'infoHanso', infoHanso ); 

end
