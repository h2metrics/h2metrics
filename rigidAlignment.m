function [dAligned, gaOpt] = rigidAlignment( dList, splineData, ...
                                             quadData, varargin) 

optTra = true;
optRot = true;
optShift = false; % Constant shifts of the parametrization
useComp = false;
maxIter = [];
display = 'off';

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
if isfield(options, 'optTra')
    optTra = options.optTra;
end
if isfield(options, 'optRot')
    optRot = options.optRot;
end
if isfield(options, 'optShift')
    optShift= options.optShift;
end
if isfield(options, 'rigidUseComp')
    useComp = options.rigidUseComp;
end
if isfield(options, 'rigidMaxIter')
    maxIter = options.rigidMaxIter;
end
if isfield(options, 'rigidDisplay')
    display = options.rigidDisplay;
end

dSpace = splineData.dSpace;
N = splineData.N;

noCurves = length(dList);
dAligned = dList;
gaOpt = cell([noCurves, 1]);

if noCurves < 2
    return
end

%% Optimization settings
options = optimoptions('fminunc');
options = optimoptions(options,'Display', display);
options = optimoptions(options,'DerivativeCheck', 'off');
% options = optimoptions(options,'PlotFcns', @optimplotfval);
options = optimoptions(options,'GradObj', 'off');
options = optimoptions(options,'Hessian', 'off');
options = optimoptions(options,'Algorithm', 'quasi-newton');

options = optimoptions(options, 'MaxFunEvals',300000);
options = optimoptions(options, 'TolFun', 1e-6);
options = optimoptions(options, 'TolX', 1e-6);
if ~isempty(maxIter) 
    options = optimoptions(options, 'maxIter', maxIter);
end

%% Optimization
F = @(coefs) rigidAlignmentDist( coefs(1:noCurves-1), ...
    coefs(noCurves:2*(noCurves-1)), ...
    reshape(coefs(2*(noCurves-1)+1:end), [noCurves-1, dSpace]), ...
    dAligned(1:noCurves), splineData, quadData, ...
    optTra, optRot, optShift );

init_coefs = zeros([(2+dSpace)*(noCurves-1), 1]);

optimal_coefs = fminunc( F, init_coefs, options );

%% Apply transformations to curve
alpha = [0; optimal_coefs(1:noCurves-1)];
beta = [0; optimal_coefs(noCurves:2*(noCurves-1))];
lambda = [0,0; reshape(optimal_coefs(2*(noCurves-1)+1:end), ...
               [noCurves-1, dSpace])];

for jj = 1:noCurves
    gaOpt{jj} = struct( 'phi', [], 'beta', [], 'v', [], 'alpha', []);
    if optTra
        dAligned{jj} = dList{jj} + ...
            ones([N, 1]) * lambda(jj,:);
        gaOpt{jj}.v = lambda(jj,:)';
    end
    if optRot
        rotation = [ cos(beta(jj)), sin(beta(jj)); ...
                    -sin(beta(jj)), cos(beta(jj)) ];
        dAligned{jj} = dAligned{jj} * rotation;
        gaOpt{jj}.beta = beta(jj);
    end
    
    if optShift
        gaOpt{jj}.alpha = alpha(jj);
        if useComp
            dAligned{jj} = curveApplyShift( dAligned{jj}, alpha(jj), ...
                                            splineData, quadData );
        else
            % We assume uniform knots and simply shift the control 
            % point sequence
            shift = round(-alpha(jj) * N / (2*pi));
            dAligned{jj} = circshift(dAligned{jj}, [-shift, 0]);
        end
    end
end

end


function D = rigidAlignmentDist( alpha, beta, lambda, ...
                                 dList, splineData, quadData, ...
                                 optTra, optRot, optShift )
% Only for dSpace=2 and periodic curves

dSpace = size(dList{1}, 2);
N = splineData.N;
nS = splineData.nS;
knotsS = splineData.knotsS;
quadPointsS = quadData.quadPointsS;
quadWeightsS = quadData.quadWeightsS;
noQuadPointsS = quadData.noQuadPointsS;

noCurves = length(dList);
dTransformed = {};

%% No transformation is applied to first curve
alpha = [0; alpha];
beta = [0; beta];
lambda = [ zeros([1, dSpace]); lambda];
                                 
%% Apply translations (lambda) and rotations (beta)
for jj = noCurves:-1:1
    if optTra
        dTransformed{jj} = dList{jj} + ones([N, 1]) * lambda(jj,:);
    end
    if optRot
        rotation = [ cos(beta(jj)), sin(beta(jj)); ...
                     -sin(beta(jj)), cos(beta(jj)) ];
        dTransformed{jj} = dTransformed{jj} * rotation;
    end
end

%% Evaluate shifted (alpha) curves at quadrature sites
cList = {};
for jj = noCurves:-1:1
    if optShift
        cList{jj} = deBoor( knotsS, nS, dTransformed{jj}, ...
                            mod(quadPointsS - alpha(jj), 2*pi), 1, ...
                            'periodic', true );
    else
        cList{jj} = quadData.B_S * dTransformed{jj};
    end
end

%% Compute distance
D = 0;
for jj = 1:noCurves-1
    for kk = jj+1:noCurves
        D = D + sqrt( sum( sum( ...
            (cList{jj}- cList{kk}).^2, 2) .* quadWeightsS));
    end
end

end
