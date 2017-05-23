%% geodesicBvp
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
%           optTra = {true, false (default)}
%           optRot = {true, false (default)}
%           optDiff = false
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
function [optE, optPath, optGa, info] = geodesicBvpParam(d0, d1, ...
    splineData, options, varargin)

%% Default parameters
optTra = true;
optRot = true;
optShift = false; % Constant shifts of the parametrization

dInitPath = [];
initGa = [];

%% Extract parameters
N = splineData.N;
Nt = splineData.Nt;
nS = splineData.nS;
nT = splineData.Nt;
dSpace = splineData.dSpace;

%% Some code for handling optional inputs
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
    ii = ii + 1;  
    end
end

%% Set options
if isfield(options, 'optTra')
    optTra = options.optTra;
end
if isfield(options, 'optRot')
    optRot = options.optRot;
end
if isfield(options, 'optShift')
    optShift= options.optShift;
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
if isempty(initGa)
    [~, gaTmp] = rigidAlignment({d0, d1}, splineData, 'options', options);
    initGa = gaTmp{2};
end

if isempty(dInitPath)
    d1Ga = curveApplyGamma(d1, initGa, splineData);
    dInitPath = linearPath(d0, d1Ga, splineData);
end

%% Setup optimization
Nphi = 1;
coeffInit = zeros([ N*(Nt-2)*dSpace + Nphi + dSpace + 2, 1]);
coeffInit(1:N*(Nt-2)*dSpace) = reshape( dInitPath(N+1:end-N, :), ...
                                        [N*(Nt-2)*dSpace, 1] );
coeffInit(end-Nphi-dSpace-2+1:end-dSpace-2) = 0; % phi
coeffInit(end-dSpace-2+1:end-2) = zeros([ dSpace, 1]); % Translation
coeffInit(end-1) = 0; % Rotation
coeffInit(end) = 0; % Shift
if optTra
    coeffInit(end-dSpace-2+1:end-2) = initGa.v;
end
if optRot
    coeffInit(end-1) = initGa.beta;
end
if optShift
    coeffInit(end) = initGa.alpha;
end

Fopt = @(coeff) energyH2Diff( ...
    [d0; reshape(coeff(1:N*(Nt-2)*dSpace), [N*(Nt-2), dSpace]); d1], ...
    coeff(end-Nphi-dSpace-2+1:end-dSpace-2), ...
    coeff(end-dSpace-2+1:end-2), coeff(end-1), coeff(end), ...
    splineData, 'optDiff', false, 'optTra', optTra, 'optRot', optRot, ...
    'optShift', optShift );

problem = struct( 'objective', Fopt, 'x0', coeffInit, ...
                  'options', minOptions, 'solver', 'fminunc' );
                  
%% Optimize
% tic
[coeffOptimal, optE, exitflag, output] = fminunc( problem );
% toc

%% Save results
% Create transformation struct
optGa = struct( 'phi', [], 'beta', [], 'v', [], 'alpha', []);
dEnd = d1;
if optShift
    optGa.alpha = coeffOptimal(end);
    dEnd = curveApplyShift(dEnd, optGa.alpha, splineData);
end
if optTra
    optGa.v = coeffOptimal(end-dSpace-2+1:end-2);
    dEnd = dEnd + ones([N, 1]) * optGa.v';
end
if optRot
    optGa.beta = coeffOptimal(end-1);
    rotation = [ cos(optGa.beta), sin(optGa.beta); ...
                 -sin(optGa.beta), cos(optGa.beta) ];
    dEnd = dEnd * rotation;
end

optPath = [ d0; ...
          reshape(coeffOptimal(1:N*(Nt-2)*dSpace), [N*(Nt-2), dSpace]); ...
          dEnd ];
      
info = struct( 'exitFlag', exitflag, ...
               'noIter', output.iterations );

end
