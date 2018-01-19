%% geodesicBvpDiff
%
% Calculates the minimal geodesic betweeen the curve d0 and the
% orbit d1 o Diff(S^1).
% 
% Input
%   d0, d1
%       Initial and final curves. Matrix of dimensions [N, dSpace].
%   splineData
%       General information about the splines used.
%
% Output
%   dPath
%       Optimal path between d0 and d1 o psi
%   psi
%       Optimal reparametrization of d1
%
function [optE, dPath, optGa, info] = geodesicBvpDiff(d0, d1, ...
    splineData, options, varargin)

% Default values
optDiff = true;
optShift = true; % Constant shifts of the parametrization -- always used,
                 % if optDiff = true;
optTra = false;
optRot = false; 

dInitPath = [];
initGa = [];

%% Extract parameters
N = splineData.N;
Nt = splineData.Nt;
nS = splineData.nS;
nT = splineData.Nt;
dSpace = splineData.dSpace;

Nphi = splineData.Nphi;
nPhi = splineData.nPhi;

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
    end
    ii = ii + 1;
end

%% Set options
if isfield(options, 'optTra')
    optTra = options.optTra;
end
if isfield(options, 'optRot')
    optRot = options.optRot;
end
options.optShift = true;
   
minOptions = optimoptions('fmincon');
% minOptions = optimoptions(minOptions,'Algorithm', ...
%                                      'trust-region-reflective');
minOptions = optimoptions(minOptions,'Algorithm', 'interior-point');

Hopt = @(coeff, lambda) energyH2DiffHessian( ...
    [d0; reshape(coeff(1:N*(Nt-2)*dSpace), [N*(Nt-2), dSpace]); d1], ...
    coeff(end-Nphi-dSpace-2+1:end-dSpace-2), ...
    coeff(end-dSpace-2+1:end-2), coeff(end-1), coeff(end), ...
    splineData, ...
    'optDiff', optDiff, 'optTra', optTra, ...
    'optRot', optRot, 'optShift', optShift );

minOptions = optimoptions(minOptions,'Hessian', 'user-supplied');
minOptions = optimoptions(minOptions,'HessFcn', Hopt);
    
minOptions = optimoptions(minOptions,'DerivativeCheck', 'off');
% minOptions = optimoptions(minOptions,'PlotFcns', @optimplotfval);
minOptions = optimoptions(minOptions,'GradObj', 'on');
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

%% Output function
fvalList = [];
    function stop = outfun(~, optimValues, ~)
        fvalList(end+1) = optimValues.fval;
        stop = false;
    end
minOptions = optimoptions(minOptions, 'OutputFcn', @outfun);

%% Generate constraints

% As phi = Id + f, the constraints encode that the control points of Id+f
% have to be increasing by at least phiEps
d_greville = aveknt(splineData.knotsPhi, nPhi+1)'; % Control points of Id

A_diff = zeros([Nphi, N*dSpace*(Nt-2)+Nphi+dSpace+2]);
for kk = 1:Nphi-1
    A_diff(kk, N*dSpace*(Nt-2) + kk) = 1;
    A_diff(kk, N*dSpace*(Nt-2) + kk + 1) = -1;
end
A_diff(Nphi, N*dSpace*(Nt-2) + Nphi) = 1;
A_diff(Nphi, N*dSpace*(Nt-2) + 1) = -1;

A_diff = sparse(A_diff);

b_diff = diff(d_greville);
b_diff = b_diff(1:end-(nPhi-1)) - splineData.phiEps;

%% Equality constraints -- average displacement of phi equals zero

Aeq = sparse(ones(1,Nphi),N*dSpace*(Nt-2)+(1:Nphi),...
    ones(Nphi,1),1,N*dSpace*(Nt-2)+Nphi+dSpace+2);
    

%% Create initial guess for path if not provided one
if isempty(initGa)
    [~, gaTmp] = rigidAlignment({d0, d1}, splineData, 'options', options);
    initGa = gaTmp{2};
    initGa.phi = zeros([ Nphi, 1 ]);
else
    if isempty(initGa.phi)
        initGa.phi = zeros([ Nphi, 1 ]);
    end
end
if isempty(initGa.alpha)
    initGa.alpha = 0;
end

if isempty(dInitPath)
    d1Ga = curveApplyGamma(d1, initGa, splineData);
    dInitPath = linearPath(d0, d1Ga, splineData);
end

%% Setup optimization
coeffInit = zeros([ N*(Nt-2)*dSpace + Nphi + dSpace + 2, 1]);
coeffInit(1:N*(Nt-2)*dSpace) = reshape( dInitPath(N+1:end-N, :), ...
                                        [N*(Nt-2)*dSpace, 1] );
coeffInit(end-Nphi-dSpace-2+1:end-dSpace-2) = initGa.phi; % phi
coeffInit(end-dSpace-2+1:end-2) = zeros([ dSpace, 1]); % Translation
coeffInit(end-1) = 0; % Rotation
coeffInit(end) = initGa.alpha; % Shift
if optTra
    coeffInit(end-dSpace-2+1:end-2) = initGa.v;
end
if optRot
    coeffInit(end-1) = initGa.beta;
end

Fopt = @(coeff) energyH2Diff( ...
    [d0; reshape(coeff(1:N*(Nt-2)*dSpace), [N*(Nt-2), dSpace]); d1], ...
    coeff(end-Nphi-dSpace-2+1:end-dSpace-2), ...
    coeff(end-dSpace-2+1:end-2), coeff(end-1), coeff(end), ...
    splineData, 'optDiff', optDiff, 'optTra', optTra, ...
    'optRot', optRot, 'optShift', optShift );

problem = struct( 'objective', Fopt, 'x0', coeffInit, ...
                  'Aineq', A_diff, 'bineq', b_diff, ...
                  'Aeq', Aeq, 'beq', 0, ...
                  'options', minOptions, 'solver', 'fmincon' );
[coeffOptimal, optE, exitflag, output] = fmincon( problem );

%% Create output
% Transformation struct
optGa = struct( 'phi', [], 'beta', [], 'v', [], 'alpha', []);

optGa.phi = coeffOptimal(end-Nphi-dSpace-2+1:end-dSpace-2);
optGa.alpha = coeffOptimal(end);
if optTra
    optGa.v = coeffOptimal(end-dSpace-2+1:end-2);
end
if optRot
    optGa.beta = coeffOptimal(end-1);
end
dEnd = curveApplyGamma(d1, optGa, splineData);

dPath = [ d0; ...
          reshape(coeffOptimal(1:N*(Nt-2)*dSpace), [N*(Nt-2), dSpace]); ...
          dEnd ];
      
fvalList = fvalList(2:end-1); % First and last entries are doubled.
      
info = struct( 'exitFlag', exitflag, ...
               'noIter', output.iterations, ...
               'fvalList', fvalList ); 

end

function [H] = energyH2DiffHessian( dPath, phi, v, beta, alpha, ...
    splineData, varargin)

[~,~,H] = energyH2Diff(dPath, phi, v, beta, alpha,...
    splineData, varargin{:} );

end
