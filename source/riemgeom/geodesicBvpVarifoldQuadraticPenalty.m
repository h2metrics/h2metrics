 %% geodesicBvpVarifoldQuadraticPenalty
%
% Minimizes the H^2 inexact boundary value problem with varifold data
% attachment term using the HANSO limited-memory BFGS method. Finds the
% geodesic between
%     d0 and d1 \circ G
% where G is a group of transformations. This can be any combination of
%   - Translations
%   - Rotations
%   - Scalings (requires scale invariant metrics)
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
%           optScal = {true, false (default)}
%           varLambda
%           hansNormTol = 1e-3 (default)
%           hansoMaxIt = 1000 (default)
%           hansoNvec = 500 (default)
%           hansoPrtLevel = {0, 1 (default), 2}
%   initPath
%       Guess for initial path. Path needs to have the same splineData as 
%       splineData Default initial path is constant path d0.
%   gaInit
%       Guess for initial transformation of d1.
%
% Output
%   optE
%       Final value of optimization routine
%   optPath
%       Optimal path between d0 and d1 o optGa
%   optGa
%       Transformation between d1 and endpoint of optPath
%   info
%       Structure containing information about the minimization process
%
function [optE, optPath, optGa, info] = geodesicBvpVarifoldQuadraticPenalty(d0, d1, ...
    splineData, options, varargin)

%% Default parameters
initPath = [];
initGa = struct('rho', [], 'beta', [], 'v', []);

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


if ~isempty(splineData.scaleInv)
    scaleInv = splineData.scaleInv;
else
    scaleInv = false;
end

if isfield( options, 'optScal' )
    optScal = options.optScal;
    if optScal == true && scaleInv == false 
       disp('optScal is only a valid option for scale invariant metrics.') 
       return
    end     
else
    optScal = false;
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

if isfield(initGa, 'rho') && ~isempty(initGa.rho)
    initRho = initGa.rho;
else
    initRho = 1;
end

coeffInit = [ reshape(initPath(splineData.N+1:end, :), [], 1); ...
              initRho;initBeta;initV ];

%% Setup HANSO
Fopt = @(coeff, pars) energyH2VarifoldQuadraticPenalty(coeff, pars, pars.splineData);

pars = struct();
pars.nvar = length(coeffInit);
pars.fgname = Fopt; %[f,df] = fgtest(x,pars)
pars.splineData = splineData;
pars.optRot = optRot;
pars.optTra = optTra;
pars.optScal = optScal;
pars.d0 = d0;
pars.dEnd = d1;
pars.lambda = options.varLambda;

if isfield(options, 'checkTurningNumber')
    pars.checkTurningNumber =  options.checkTurningNumber;
else
    pars.checkTurningNumber = 0;
end

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
if isfield(options, 'hansoCpuMax')
    optionsHANSO.cpumax = options.hansoCpuMax;
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

% Ask Hanso to continue even when the line search fails.
% This is used in conjunction with energy=inf in case of changing turning
% number.
optionsHANSO.quitLSfail = 0;

%% Call HANSO
[optCoeff, optE, infoHanso] = hanso(pars, optionsHANSO);

%% Create output
% Transformation struct
optGa = struct( 'rho',[],'beta', [], 'v', [] );

if optTra
    optGa.v = optCoeff(end-dSpace+1:end);
end
if optRot
    optGa.beta = optCoeff(end-dSpace);
end

if optScal
    optGa.rho = optCoeff(end-dSpace-1);
end

optPath = [ d0; reshape( optCoeff(1:end-dSpace-2), [], 2) ];
dPathEnd = optPath(end-splineData.N+1:end,:);
dEnd = curveApplyGamma( d1, optGa, splineData);

% Construct derivative path similar to optPath
grad = infoHanso.grad;
dOptPath = [zeros(size(d0)); reshape( grad(1:end-dSpace-2), [], 2) ];
infoHanso.grad = dOptPath;

distVarSqrd = varifoldDistanceSquared(dPathEnd, dEnd, splineData);
    
optE = optE - options.varLambda*distVarSqrd;

info = struct( 'infoHanso', infoHanso, 'optDistVar',distVarSqrd ); 

end
