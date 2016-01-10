function [dMean, info] = karcherMeanManopt(dList, ...
    splineData, quadData, quadDataTensor, varargin)

options = [];
dInit = [];

% Some code for handling optional inputs
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'options'
                ii = ii + 1;
                options = varargin{ii};
            case 'meaninit'
                ii = ii + 1;
                dInit = varargin{ii};
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1;  
    end
end

N = splineData.N;
dSpace = splineData.dSpace;
noCurves = length(dList);

splineData.stepsT = splineData.stepsT;

curveManifold = constructCurveManifold(splineData, quadData);

dLast = [];
dLastPathList = [];
gaLastList = [];

function stats = statsfun(problem, dIter, stats, store)
    if ~isfield(store, 'ready')
        [~, store] = cost(dIter, store);
    end
    
    dLast = dIter;
    dLastPathList = store.dPathList;
    gaLastList = store.gaList;
    
    stats.dIter = dIter;
    stats.dPathList = store.dPathList;
    stats.gaList = store.gaList;
    stats.enList = store.enList;
end

function [sumE, store] = cost(dIter, store)
    if ~isfield(store, 'ready')
        dPathList = {};
        gaList = {};
        enList = [];
        for jj = noCurves:-1:1
            disp(jj);
            if ~isempty(dLast)
                dInitPath = dLastPathList{jj} + ...
                    linearPath( dIter - dLast, zeros([N, dSpace]), ...
                                splineData );
                
                [optE, dPath, optGa] = geodesicBvpDiff(dIter, dList{jj},...
                    splineData, quadData, quadDataTensor, ...
                    'options', options, 'initPath', dInitPath, ...
                    'initGa', gaLastList{jj});
            else
                [optE, dPath, optGa] = ...
                    geodesicBvpDiff(dIter, dList{jj}, splineData, ...
                        quadData, quadDataTensor, 'options', options);
            end
            dPathList{jj} = dPath;
            gaList{jj} = optGa;
            enList(jj) = optE;
        end
        store.dPathList = dPathList;
        store.gaList = gaList;
        store.enList = enList;
        store.ready = true;
    end
    
    sumE = sum(store.enList) / noCurves;
end

function [g, store] = grad(dIter, store)
    if ~isfield(store, 'ready')
        [~, store] = cost(dIter, store);
    end
    
    g = zeros([N, dSpace]);
    for jj = noCurves:-1:1
        g = g + pathVelocity(store.dPathList{jj}, 0, splineData);
    end
    g = -2 * g / noCurves;
end

% Setup the problem structure
problem.M = curveManifold;
problem.cost = @cost;
problem.grad = @grad;

% Initial guess
if isempty(dInit)
    dInit = dList{1};
end

% Solver options
karcherOptions = [];
if isfield(options, 'karcherTolGradNorm')
    karcherOptions.tolgradnorm = options.karcherTolGradNorm;
end
if isfield(options, 'karcherMaxIter')
    karcherOptions.maxiter = options.karcherMaxIter;
end

karcherOptions.statsfun = @statsfun;
karcherOptions.verbosity = 3;

[dMean, ~, info] = conjugategradient(problem, dInit, karcherOptions);
% [dMean, E, info] = steepestdescent(problem, dInit, opts);
    
end

function M = constructCurveManifold( splineData, quadData )
    
    M.splineData = splineData;
    M.quadData = quadData;
    
    M.name = @() sprintf('Curve manifold with N=%d, n=%d', ...
        splineData.N, splineData.nS);

    M.dim = @() splineData.N * splineData.dSpace;
    
    M.inner = @(c, v1, v2) ...
        curveRiemH2InnerProd(c, v1, v2, splineData, quadData);
    
    M.norm = @(c, v) sqrt(M.inner(c, v, v));
    
    % Expensive to compute
    % M.dist = @(x, y) real(acos(x(:).'*y(:)));
    
    % IMPROVE THIS
    M.typicaldist = @() 1;
    
    M.proj = @(c, v) v;
    
    M.tangent = M.proj;
	
    % We can implement this, if needed later
	% M.egrad2rgrad = M.proj;
	
    % Hopefully not needed
	% M.ehess2rhess = @ehess2rhess;
	% function rhess = ehess2rhess(x, egrad, ehess, u)
    %     rhess = M.proj(x, ehess) - (x(:)'*egrad(:))*u;
	% end
    
    M.exp = @(c, v, t) exponential(c, v, t, splineData, quadData);
    
    M.retr = @(c, v, t) M.exp(c, v, t);

    % M.log = @logarithm;
    % function v = logarithm(x1, x2)
    %     v = M.proj(x1, x2 - x1);
    %     di = M.dist(x1, x2);
	%     nv = norm(v, 'fro');
	%     v = v * (di / nv);
    % end
    
    M.hash = @(x) ['z' hashmd5(x(:))];
    
    % M.rand = @() random(n, m);
    
    % M.randvec = @(x) randomvec(n, m, x);
    
    M.lincomb = @lincomb;
    
    M.zerovec = @(v) zeros([splineData.N, splineData.dSpace]);
    
    M.transp = @(c1, c2, v) v;
    
    % M.pairmean = @pairmean;
    % function y = pairmean(x1, x2)
    %     y = x1+x2;
    %     y = y / norm(y, 'fro');
    % end

    M.vec = @(x, u_mat) ...
        reshape(u_mat, [splineData.N * splineData.dSpace, 1]);
    M.mat = @(x, u_vec) ...
        reshape(u_vec, [splineData.N, splineData.dSpace]);
    M.vecmatareisometries = @() false;

end

% Exponential on the sphere
function d = exponential(c, v, t, splineData, quadData)
    d = geodesicForward(c, c + t / splineData.stepsT * v, ...
        splineData.stepsT, splineData, quadData, 'endpoint');
end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2)

    if nargin == 3
        d = a1*d1;
    elseif nargin == 5
        d = a1*d1 + a2*d2;
    else
        error('Bad use of lincomb.');
    end

end
