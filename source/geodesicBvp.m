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

initPath = [];
initGa = [];

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

%% Multigrid optimization and recursive call
if isfield(options, 'useMultigrid') && options.useMultigrid
    % Create roughSplineData and roughOptions if not provided
    if isfield(options, 'mgSplineData')
        roughSplineData = options.mgSplineData;
    else
        roughSplineData = splineData;
        roughSplineData.Nt = min(3, splineData.Nt);
        roughSplineData = constructKnots(roughSplineData);
        roughSplineData = setupQuadData(roughSplineData);
    end
    if isfield(options, 'mgOptions')
        roughOptions = options.mgOptions;
    else
        roughOptions = options;
        roughOptions.useMultigrid = false;
        roughOptions.mgSplineData = [];
    end
    
    % Resample curves and initial data
    d0R = curveSpline2Spline(d0, splineData, roughSplineData);
    d1R = curveSpline2Spline(d1, splineData, roughSplineData);
    if ~isempty(initPath)
        initPathR = pathSpline2Spline(initPath, splineData, ...
                                      roughSplineData);
    else
        initPathR = [];
    end
    initGaR = initGa;
    
    % Call geodesicBvp
    [~, optPathR, optGaR, ~] = geodesicBvp(d0R, d1R, roughSplineData, ...
        roughOptions, 'initPath', initPathR, 'initGa', initGaR);
    
    % Set varargin for function calls below
    initPath = pathSpline2Spline(optPathR, roughSplineData, splineData);
    initGa = optGaR;
    
    varargin = {'initPath', initPath, 'initGa', initGa};                 
end

%% Decision tree

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
