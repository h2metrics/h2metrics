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

options = [];

% Some code for handling optional inputs
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'options'
                ii = ii + 1;
                options = varargin{ii};
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
   
if optDiff
    minOptions = optimoptions('fmincon');
    minOptions = optimoptions(minOptions,'Algorithm', 'interior-point');
else
    minOptions = optimoptions('fminunc');
    minOptions = optimoptions(minOptions,'Algorithm', 'quasi-newton');
end
minOptions = optimoptions(minOptions,'DerivativeCheck', 'off');
% minOptions = optimoptions(minOptions,'PlotFcns', @optimplotfval);
minOptions = optimoptions(minOptions,'GradObj', 'off');
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

%% Extract parameters
N = splineData.N;
Nt = splineData.Nt;
nS = splineData.nS;
nT = splineData.Nt;
Nphi = splineData.Nphi;
nPhi = splineData.nPhi;
dSpace = splineData.dSpace;

%% Generate constraints
if optDiff
    % As phi = Id + f, the constraints encode that the control points of Id+f
    % have to be increasing by at least phiEps
    d_greville = aveknt(splineData.knotsPhi, nPhi+1)'; % Control points of Id

    A_diff = zeros([Nphi+nPhi-1, N*dSpace*(Nt-2)+Nphi+dSpace+1]);
    for kk = 1:Nphi-1
        A_diff(kk, N*dSpace*(Nt-2) + kk) = 1;
        A_diff(kk, N*dSpace*(Nt-2) + kk + 1) = -1;
    end
    A(Nphi, N*dSpace*(Nt-2) + Nphi) = 1;
    A(Nphi, N*dSpace*(Nt-2) + 1) = -1;
    for kk = 1:nPhi-1 % Because of periodicity, the first control points are
                      % repeated at the end
        A_diff(Nphi + kk, N*dSpace*(Nt-2) + kk) = 1;
        A_diff(Nphi + kk, N*dSpace*(Nt-2) + kk + 1) = -1;
    end
    A_diff = sparse(A_diff);
    b_diff = diff(d_greville) + splineData.phiEps;

    % phi1_nonper = [ zeros([N*dSpace*(Nt-2), 1]); phi1 ];
    % disp(A_diff * phi1_nonper < b_diff);
else
    Nphi = 1; % Simpler than setting it 0.
end
% if ~optDiff
%     A_diff = sparse(zeros([1, N*dSpace*(Nt-2)+Nphi+dSpace+1]));
%     b_diff = 0;
% end

%% Setup optimization

dLinear = linearPath(d0, d1, splineData);
coeffInit = zeros([ N*(Nt-2)*dSpace + Nphi + dSpace + 1, 1]);
coeffInit(1:N*(Nt-2)*dSpace) = reshape( dLinear(N+1:end-N, :), ...
                                        [N*(Nt-2)*dSpace, 1] );
coeffInit(end-Nphi-dSpace-1+1:end-dSpace-1) = zeros([ Nphi, 1]); % phi
coeffInit(end-dSpace-1+1:end-1) = zeros([ dSpace, 1]); % Translation
coeffInit(end) = 0; % Rotation

Fopt = @(coeff) energyH2Diff( ...
    [d0; reshape(coeff(1:N*(Nt-2)*dSpace), [N*(Nt-2), dSpace]); d1], ...
    coeff(end-Nphi-dSpace-1+1:end-dSpace-1), ...
    coeff(end-dSpace-1+1:end-1), coeff(end), ...
    splineData, quadData, quadDataTensor, ...
    'optDiff', optDiff, 'optTra', optTra, 'optRot', optRot );

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
optGa = struct( 'phi', [], 'beta', [], 'v', [] );
dEnd = d1;
if optDiff
    optGa.phi = coeffOptimal(end-Nphi-dSpace-1+1:end-dSpace-1);
    dEnd = curveComposeDiff(dEnd, optGa.phi, splineData, quadData);
end
if optTra
    optGa.v = coeffOptimal(end-dSpace-1+1:end-1);
    dEnd = dEnd + ones([N, 1]) * optGa.v';
end
if optRot
    optGa.beta = coeffOptimal(end);
    rotation = [ cos(optGa.beta), -sin(optGa.beta); ...
                 sin(optGa.beta),  cos(optGa.beta) ];
    dEnd = dEnd * rotation;
end

dPath = [ d0; ...
          reshape(coeffOptimal(1:N*(Nt-2)*dSpace), [N*(Nt-2), dSpace]); ...
          dEnd ];
      
exitFlag = exitflag;
noIter = output.iterations;

end
