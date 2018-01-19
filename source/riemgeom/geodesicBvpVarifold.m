%% geodesicBvpVarifold
%
% Minimizes the H^2 inexact boundary value problem with varifold data
% attachment term using the HANSO limited-memory BFGS method. Finds the
% geodesic between
%     d0 and d1 \circ G
% where G is a group of transformations. This can be any combination of
%   - Translations
%   - Rotations
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
%       Final value of geodesic energy (does not include varifold terms)
%   optPath
%       Optimal path between d0 and d1 o optGa
%   optGa
%       Transformation between d0 and endpoint of optPath
%   optDistVar
%       Varifold distance between end point of optimal path and target
%       curve
%   info
%       Structure containing information about the minimization process
%
% TODO: Add optimization parameters for Aug. Lag. to splineData or 
% somewhere appropriate:
% Change Update rule for etaK? (Now *5)
%
%
function [optE, optPath, optGa, info] = geodesicBvpVarifold(d0, d1, ...
    splineData, options, varargin)

%% Default parameters
initPath = [];
initGa = struct('beta', [], 'v', [], 'rho',[]);

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
              initRho;...
              initBeta; ...
              initV ];

%% Setup HANSO
Fopt = @(coeff, pars) energyH2Varifold(coeff, pars, pars.splineData);

pars = struct();
pars.nvar = length(coeffInit);
pars.fgname = Fopt; %[f,df] = fgtest(x,pars)
pars.splineData = splineData;
pars.optRot = optRot;
pars.optTra = optTra;
pars.optScal = optScal;
pars.d0 = d0;
pars.dEnd = d1;


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

% %% Call HANSO
% [optCoeff, optE, infoHanso] = hanso(pars, optionsHANSO);

%% Run Augmented Lagrangian 
% Setup default parameters. TODO: Add to splineData.varData ? Somewhere
% appropriate
maxIter = 50; %Max no. iterations for Aug. Lag.
% eps_match = 1e-2; %Soft constraint for similarity term
% tau_final = 1e-3; %Norm gradient criteria

% Get values from splineData instead
eps_match =  options.eps_match;
tau_final =  options.tau_final;

lambdaK = zeros(1,maxIter); %Initial value of Lagrangian 
etaK = zeros(1,maxIter); etaK(1) = 1; 
tauK = zeros(1,maxIter); tauK(1) = 1e-1;

optEK = zeros(1,maxIter);
distVarK = zeros(1,maxIter);

for k = 1:maxIter
    % Call Hanso
    pars.lambda = lambdaK(k);
    pars.eta = etaK(k);
    optionsHANSO.normtol = tauK(k);
        
    [optCoeff, optL, infoHanso] = hanso(pars, optionsHANSO);
    
    % Extract end curve and compute varifold distance
    v = zeros(2,1);
    beta = 0;
    rho = 1;
    if optTra
        v = optCoeff(end-dSpace+1:end);
    end
    if optRot
        beta = optCoeff(end-dSpace);
    end
    if optScal
        rho = optCoeff(end-dSpace-1);
    end
    optPath = [ d0; reshape( optCoeff(1:end-dSpace-2), [], 2) ];
    dPathEnd = optPath(end-splineData.N+1:end,:);
    rotation = [ cos(beta), sin(beta); ...
        -sin(beta), cos(beta) ];
    dEnd = d1 + ones([splineData.N, 1]) * v(:)';
    dEnd = dEnd * rotation;
    dEnd = dEnd * rho;
    
    
    distVarSqrd = varifoldDistanceSquared(dPathEnd, dEnd, splineData);
    
    % Test convergence
    if distVarSqrd > eps_match %Constraint violated
        % Update lambda
        lambdaK(k+1) = lambdaK(k) - etaK(k)*distVarSqrd;
        etaK(k+1) = 5*etaK(k);
    else % Constraints are ok
        if tauK(k) <= tau_final % Converged
            % Save results
            optEK(k) = optL + lambdaK(k)*distVarSqrd - etaK(k)/2*distVarSqrd^2;
            distVarK(k) = distVarSqrd;
            break
        else 
            lambdaK(k+1) =  lambdaK(k);
            etaK(k+1) = etaK(k);
        end
    end
    if tauK(k) > tau_final
        tauK(k+1) = tauK(k)/10;
    else
        tauK(k+1) = tauK(k);
    end
    
    % Update initial point
    optionsHANSO.x0 = optCoeff;
    
    % Save results
    optEK(k) = optL + lambdaK(k)*distVarSqrd - etaK(k)/2*distVarSqrd^2;
    distVarK(k) = distVarSqrd;
end

if k == maxIter
    disp('Warning: Max number of iterations in augmented Lagrangian reached')
end

%% Create output
% Transformation struct
%TODO: the break completely puts us out of the function!
optGa = struct('rho',[], 'beta', [], 'v', [] );
% if optTra
%     optGa.v = optCoeff(end-dSpace+1:end);
% end
% if optRot
%     optGa.beta = optCoeff(end-dSpace);
% end
% if optScal
%     optGa.rho = optCoeff(end-dSpace-1);
% end

%Update the gamma structure
optGa.v = optCoeff(end-dSpace+1:end);
optGa.beta = optCoeff(end-dSpace);
optGa.rho = optCoeff(end-dSpace-1);

optPath = [ d0; reshape( optCoeff(1:end-dSpace-2), [], 2) ];
           
info = struct( 'infoHanso', infoHanso); 

% Compute the final squared varifold distance
% Extract end points, and compute varifold distance.
rotation = [ cos(optGa.beta), sin(optGa.beta); ...
    -sin(optGa.beta), cos(optGa.beta) ];
dEnd = d1 + ones([splineData.N, 1]) * optGa.v(:)';
dEnd = dEnd * rotation;
dEnd = dEnd * optGa.rho;

dPathEnd = optPath(end-splineData.N+1:end,:);
distVarSqrd = varifoldDistanceSquared(dPathEnd, dEnd, splineData);
info.optDistVar = distVarSqrd;
info.dEnd = dEnd;
info.dPathEnd = dPathEnd;
% Compute the Riemannian Energy of the optimal path
optE = optL + lambdaK(k)*distVarSqrd - etaK(k)/2*distVarSqrd^2;

info.lambdaK = lambdaK(1:k);
info.etaK = etaK(1:k);
info.tauK = tauK(1:k);

info.optEK = optEK(1:k);
info.distVarK = distVarK(1:k);

end
