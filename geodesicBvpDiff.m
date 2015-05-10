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
function [dPath, psi] = geodesicBvpDiff(d0, d1, ...
    splineData, quadData, quadDataTensor, varargin)

options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'DerivativeCheck', 'off');
options = optimoptions(options,'PlotFcns', @optimplotfval);
options = optimoptions(options,'GradObj', 'off');
options = optimoptions(options,'Hessian', 'off');
options = optimoptions(options,'Algorithm', 'interior-point');
options = optimoptions(options,'UseParallel',false);
options = optimoptions(options,'MaxFunEvals',300000);
options = optimoptions(options,'TolFun', 1e-3);
options = optimoptions(options,'TolX', 1e-3);

% Some code for handling optional inputs
ii = 1;
while ii <= length(varargin)-1
    if (isa(varargin{ii},'char') && isa(varargin{ii+1},'char'))
         switch (lower(varargin{ii}))
             case 'display' %Don't return gradient terms from end curves
                 options = optimoptions(options,'Display',lower(varargin{ii+1}));
             case 'plotfval' %value = 'true'/'false'
                 if strcmpi(varargin{ii+1},'true')
                    options = optimoptions(options,'PlotFcns', @optimplotfval); 
                 end
             case 'algorithm'
                 options = optimoptions(options,'Algorithm',lower(varargin{ii+1}));
             case 'gradobj'
                 options = optimoptions(options,'GradObj', lower(varargin{ii+1}));
         end
    elseif (isa(varargin{ii},'char') && isa(varargin{ii+1},'numeric'))
        switch (lower(varargin{ii}))
            case 'tolfun'
                options = optimoptions(options,'TolFun', varargin{ii+1});
            case 'tolx'
                options = optimoptions(options,'TolX', varargin{ii+1});
            case 'init'
                d_init = varargin{ii+1};
            case 'MaxFunEvals'
                options = optimoptions(options, 'MaxFunEvals',varargin{ii+1});
        end
    end
    ii = ii + 1;  
end

%% Extract parameters
N = splineData.N;
Nt = splineData.Nt;
Nphi = splineData.Nphi;
nS = splineData.nS;
nT = splineData.Nt;
nPhi = splineData.nPhi;
dSpace = splineData.dSpace;

%% Generate constraints
% As phi = Id + f, the constraints encode that the control points of Id+f
% have to be increasing by at least phiEps
d_greville = aveknt(splineData.knotsPhi, nPhi+1)'; % Control points of Id

A_diff = zeros([Nphi+nPhi-1, N*dSpace*(Nt-2)+Nphi]);
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

%% Setup optimization

dLinear = linearPath(d0, d1, splineData);
coeffInit = zeros([ N*(Nt-2)*dSpace + Nphi, 1]);
coeffInit(1:N*(Nt-2)*dSpace) = reshape( dLinear(N+1:end-N, :), ...
                                        [N*(Nt-2)*dSpace, 1] );
coeffInit(end-Nphi+1:end) = zeros([ Nphi, 1]); % Identity diffeomorphism

Fopt = @(coeff) energyH2Diff( ...
    [d0; reshape(coeff(1:end-Nphi), [N*(Nt-2), dSpace]); d1], ...
    coeff(end-Nphi+1:end), ...
    splineData, quadData, quadDataTensor );


problem = struct( 'objective', Fopt, 'x0', coeffInit, ...
                  'Aineq', A_diff, 'bineq', b_diff, ...
                  'options', options, 'solver', 'fmincon' );

tic
coeffOptimal = fmincon( problem );
% [coeffOptimal, EOptimal, exitflag, output] = fmincon( problem );
toc

psi = coeffOptimal(end-Nphi+1:end);
dPath = [d0; reshape(coeffOptimal(1:end-Nphi), [N*(Nt-2), dSpace])];
dPath = [ dPath; composeCurveDiff(d1, psi, splineData, quadData) ];

end

