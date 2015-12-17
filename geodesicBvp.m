%% geodesicBvp
%
% Calculates the minimal geodesic betweeen 
%     d0 and d1 o G
% where d0, d1 are curves and G is a group of transformations. At the
% moment we support G to be any combination of
%   - Translations
%   - Rotations
%   - Reparametrizations (full diffeomorphism group)
%   - Rigid shifts of the parametrization (S^1 as a subgroup of Diff(S^1))
% 
% Input
%   d0, d1
%       Initial and final curves. Matrix of dimensions [N, dSpace].
%   splineData
%       General information about the splines used.
%   quadData, quadDataTensor
%       Precomputed spline collocation matrices at quadrature points.
%
% Optional parameters
%   options
%       Struct containing optimization options. Will be passed on to
%       subroutines. Contains the following fields:
%           useAmpl = {true, false (default)}
%           optTra = {true, false (default)}
%           optRot = {true, false (default)}
%           optDiff = {true, false (default)}
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
function [optE, optPath, optGa, info] = geodesicBvp(d0, d1, ...
    splineData, quadData, quadDataTensor, varargin)

% Some code for handling optional inputs
ii = 1;
while ii <= length(varargin)
    if strcmpi(varargin{ii}, 'options')
        ii = ii + 1;
        options = varargin{ii};
    end
    ii = ii + 1;  
end

% Enforce default options
if ~isfield(options, 'useAmpl')
    options.useAmpl = false;
end

% Decision tree
% if useAmpl
%   call geodesicBvpAmpl
% elseif optDiff
%   call geodesicBvpDiff
% else
%   call geodesicBvpParam

if isfield(options, 'useAmpl') && options.useAmpl
    error('Ampl version not implemented yet...');
elseif isfield(options, 'optDiff') && options.optDiff
    [optE, optPath, optGa, info] = geodesicBvpDiff(d0, d1, ...
        splineData, quadData, quadDataTensor, varargin{:});
else
    [optE, optPath, optGa, info] = geodesicBvpParam(d0, d1, ...
        splineData, quadData, quadDataTensor, varargin{:});
end

end
