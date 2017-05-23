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
%       Initial and final curves. Matrices of dimension [N, dSpace].
%   splineData
%       General information about the splines used.
%   options
%       Struct containing optimization options. Will be passed on to
%       subroutines. Contains the following fields:
%           optTra = {true, false (default)}
%           optRot = {true, false (default)}
%           optDiff = {true, false (default)}
%           optShift = {true, false (default)}
%           maxIter = integer ([] for default value)
%           display = string
%               'off' for no output
%               '' for default of optimization routine
%               Any other string will be passed on
%
% Optional parameters
%   initPath
%       Guess for initial path.
%   gaInit
%       Guess for initial transformation of d1.
%
% Output
%   optE
%       Final value of optimization routine; usually corresponds to the
%       energy of the path.
%   optPath
%       Optimal path between d0 and d1 o G
%   optGa
%       Transformation between d1 and endpoint of optPath
%   info
%       Structure containing information about the minimization process
%
function [optE, optPath, optGa, info] = geodesicBvp(d0, d1, ...
    splineData, options, varargin)

% Decision tree
% if optDiff && useVarifold
%   call geodesicBvpVarifold
% elseif optDiff
%   call geodesicBvpDiff
% elseif noGroupsAtAll
%   call geodesicBvpNoGroups
% else
%   call geodesicBvpParam

useVarifold = isfield(options, 'useVarifold') && options.useVarifold;
doDiff = isfield(options, 'optDiff') && options.optDiff;
doTra = isfield(options, 'optTra') && options.optTra;
doRot = isfield(options, 'optRot') && options.optRot;
doShift = isfield(options, 'optShift') && options.optShift;

if doDiff && useVarifold
    [optE, optPath, optGa, info] = geodesicBvpVarifold(d0, d1, ...
        splineData, options, varargin{:});
elseif doDiff
    [optE, optPath, optGa, info] = geodesicBvpDiff(d0, d1, ...
        splineData, options, varargin{:});
elseif ~doDiff && ~doTra && ~doRot && ~doShift
    [optE, optPath, optGa, info] = geodesicBvpNoGroups(d0, d1, ...
        splineData, options, varargin{:});
else
    [optE, optPath, optGa, info] = geodesicBvpParam(d0, d1, ...
        splineData, options, varargin{:});
end

end
