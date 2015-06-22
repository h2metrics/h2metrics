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
%   quadData, quadDataTensor
%       Precomputed spline collocation matrices at quadrature points.
%
% Output
%   dPath
%       Optimal path between d0 and d1 o psi
%   psi
%       Optimal reparametrization of d1
%
function [optE, dPath, optGa, exitFlag, noIter] = geodesicBvpDiff(d0, d1, ...
    splineData, quadData, quadDataTensor, varargin)

optDiff = true;
optTra = true;
optRot = true;
optShift = false; % Constant shifts of the parametrization

options = [];

dInitPath = [];
gaInit = [];

%% Extract parameters
N = splineData.N;
Nt = splineData.Nt;
nS = splineData.nS;
nT = splineData.Nt;
dSpace = splineData.dSpace;
if optDiff
    Nphi = splineData.Nphi;
    nPhi = splineData.nPhi;
end

% Some code for handling optional inputs
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
            case 'gainit'
                ii = ii + 1;
                gaInit = varargin{ii};
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1;  
    end
end

% Set options
if isfield(options, 'optDiff')
    optDiff = options.optDiff;
end
if isfield(options, 'optTra')
    optTra = options.optTra;
end
if isfield(options, 'optRot')
    optRot = options.optRot;
end
if isfield(options, 'optShift')
    optShift= options.optShift;
end
   
if optDiff
    minOptions = optimoptions('fmincon');
%     minOptions = optimoptions(minOptions,'Algorithm', 'trust-region-reflective');
    minOptions = optimoptions(minOptions,'Algorithm', 'interior-point');
    
    Hopt = @(coeff,lambda) energyH2DiffHessian( ...
    [d0; reshape(coeff(1:N*(Nt-2)*dSpace), [N*(Nt-2), dSpace]); d1], ...
    coeff(end-Nphi-dSpace-2+1:end-dSpace-2), ...
    coeff(end-dSpace-2+1:end-2), coeff(end-1), coeff(end), ...
    splineData, quadData, quadDataTensor, ...
    'optDiff', optDiff, 'optTra', optTra, 'optRot', optRot, ...
    'optShift', optShift );

    minOptions = optimoptions(minOptions,'Hessian', 'user-supplied',...
        'HessFcn',Hopt);
else
    minOptions = optimoptions('fminunc');
    minOptions = optimoptions(minOptions,'Algorithm', 'trust-region');
    minOptions = optimoptions(minOptions,'Hessian', 'on');
end
minOptions = optimoptions(minOptions,'DerivativeCheck', 'off');
% minOptions = optimoptions(minOptions,'PlotFcns', @optimplotfval);
minOptions = optimoptions(minOptions,'GradObj', 'on');
% minOptions = optimoptions(minOptions,'Hessian', 'off');
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

%% Extract parameters
N = splineData.N;
Nt = splineData.Nt;
nS = splineData.nS;
nT = splineData.Nt;
dSpace = splineData.dSpace;
if optDiff
    Nphi = splineData.Nphi;
    nPhi = splineData.nPhi;
end

%% Generate constraints
if optDiff
    % As phi = Id + f, the constraints encode that the control points of Id+f
    % have to be increasing by at least phiEps
    d_greville = aveknt(splineData.knotsPhi, nPhi+1)'; % Control points of Id

%     A_diff = zeros([Nphi+(nPhi-1), N*dSpace*(Nt-2)+Nphi+dSpace+2]);
    A_diff = zeros([Nphi, N*dSpace*(Nt-2)+Nphi+dSpace+2]);
    for kk = 1:Nphi-1
        A_diff(kk, N*dSpace*(Nt-2) + kk) = 1;
        A_diff(kk, N*dSpace*(Nt-2) + kk + 1) = -1;
    end
    A_diff(Nphi, N*dSpace*(Nt-2) + Nphi) = 1;
    A_diff(Nphi, N*dSpace*(Nt-2) + 1) = -1;
%     for kk = 1:nPhi-1 % Because of periodicity, the first control points are
%                       % repeated at the end
%         A_diff(Nphi + kk, N*dSpace*(Nt-2) + kk) = 1;
%         A_diff(Nphi + kk, N*dSpace*(Nt-2) + kk + 1) = -1;
%     end
    A_diff = sparse(A_diff);
    
    b_diff = diff(d_greville);
    b_diff = b_diff(1:end-(nPhi-1)) + splineData.phiEps;

    % phi1_nonper = [ zeros([N*dSpace*(Nt-2), 1]); phi1 ];
    % disp(A_diff * phi1_nonper < b_diff);
else
    Nphi = 1; % Simpler than setting to 0
end

% if ~optDiff
%     A_diff = sparse(zeros([1, N*dSpace*(Nt-2)+Nphi+dSpace+1]));
%     b_diff = 0;
% end

%% Create initial guess for path if not provided one
if isempty(gaInit)
    [~, gaTmp] = rigidAlignment( {d0, d1}, splineData, quadData, ...
                                 'options', options );
    gaInit = gaTmp{2};
end

if isempty(dInitPath)
    d1Ga = curveApplyGamma(d1, gaInit, splineData, quadData);
    dInitPath = linearPath(d0, d1Ga, splineData);
end

%% Setup optimization
coeffInit = zeros([ N*(Nt-2)*dSpace + Nphi + dSpace + 2, 1]);
coeffInit(1:N*(Nt-2)*dSpace) = reshape( dInitPath(N+1:end-N, :), ...
                                        [N*(Nt-2)*dSpace, 1] );
coeffInit(end-Nphi-dSpace-2+1:end-dSpace-2) = zeros([ Nphi, 1]); % phi
coeffInit(end-dSpace-2+1:end-2) = zeros([ dSpace, 1]); % Translation
coeffInit(end-1) = 0; % Rotation
coeffInit(end) = 0; % Shift
if optTra
    coeffInit(end-dSpace-2+1:end-2) = gaInit.v;
end
if optRot
    coeffInit(end-1) = gaInit.beta;
end
if optShift
    coeffInit(end) = gaInit.alpha;
end

Fopt = @(coeff) energyH2Diff( ...
    [d0; reshape(coeff(1:N*(Nt-2)*dSpace), [N*(Nt-2), dSpace]); d1], ...
    coeff(end-Nphi-dSpace-2+1:end-dSpace-2), ...
    coeff(end-dSpace-2+1:end-2), coeff(end-1), coeff(end), ...
    splineData, quadData, quadDataTensor, ...
    'optDiff', optDiff, 'optTra', optTra, 'optRot', optRot, ...
    'optShift', optShift );

if optDiff
    problem = struct( 'objective', Fopt, 'x0', coeffInit, ...
                      'Aineq', A_diff, 'bineq', b_diff, ...
                      'options', minOptions, 'solver', 'fmincon' );
    % tic
    [coeffOptimal, optE, exitflag, output] = fmincon( problem );
    % toc
else
    problem = struct( 'objective', Fopt, 'x0', coeffInit, ...
                      'options', minOptions, 'solver', 'fminunc' );
    % tic
    [coeffOptimal, optE, exitflag, output] = fminunc( problem );
    % toc
end

% Create transformation struct
optGa = struct( 'phi', [], 'beta', [], 'v', [], 'alpha', []);
dEnd = d1;
if optDiff && optShift
    optGa.phi = coeffOptimal(end-Nphi-dSpace-2+1:end-dSpace-2);
    optGa.alpha = coeffOptimal(end);
    % First diffeo, then shift.
    dEnd = curveComposeDiff( dEnd, optGa.phi - optGa.alpha, ...
                             splineData, quadData ); 
elseif optShift
    optGa.alpha = coeffOptimal(end);
    dEnd = curveApplyShift(dEnd, optGa.alpha, splineData, quadData);
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


dPath = [ d0; ...
          reshape(coeffOptimal(1:N*(Nt-2)*dSpace), [N*(Nt-2), dSpace]); ...
          dEnd ];
      
exitFlag = exitflag;
noIter = output.iterations;

end

function [H] = energyH2DiffHessian( dPath, phi, v, beta, alpha, ...
    splineData, quadData, quadDataTensor, varargin)

[~,~,H] = energyH2Diff(dPath, phi, v, beta, alpha,...
    splineData, quadData, quadDataTensor, varargin{:} );

end
