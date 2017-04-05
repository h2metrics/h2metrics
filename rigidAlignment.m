%% rigidAlignment
%
% Computes the rigidAlignment between a set of curves
%
% NOTE: Code works only for dSpace=2 and periodic curves. For nonperiodic
% curves shifts don't work and in more than two dimensions rotations have
% to be rethought.
%
% Input
%   dList
%       List of curves
%   splineData
%       General information about the splines used.
%
% Optional inputs
%   lineStyle = 'k-' (default)
%       lineStyle parameter to be passed to plot.
%
% Output
%   dAligned
%       List of aligned curves
%   gaOpt
%       List of gamma structures to attain alignment
%
% Notes
%   The order of precedence for the constants are as follows
%     -) Optional parameter 'a'
%     -) splineData.a
%     -) a = [1 0 1]
%
function [dAligned, gaOpt] = rigidAlignment( dList, splineData, varargin) 

optTra = true;
optRot = true;
optShift = false; % Constant shifts of the parametrization
useComp = false;
maxIter = [];
display = 'off';
globalRot = false;
a = [1 0 1];
if ~isempty(splineData.a)
    a = splineData.a;
end

options = struct;

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
if isfield(options, 'rigidGlobalRot')
    globalRot = options.rigidGlobalRot;
end
if isfield(options, 'rigidA')
    a = options.rigidA;
end
dSpace = splineData.dSpace;
N = splineData.N;

noCurves = length(dList);
dAligned = dList;
gaOpt = cell([noCurves, 1]);

if noCurves < 2
    return
end

if noCurves==2 && optShift && ~useComp && ~optTra && ~optRot
    [dAligned, gaOpt] = shiftOnlyTwoCurves(dList, splineData, a);
    return
end
    

%% Optimization settings
optionsOpt = optimoptions('fminunc');
optionsOpt = optimoptions(optionsOpt,'Display', display);
optionsOpt = optimoptions(optionsOpt,'DerivativeCheck', 'off');
% optionsOpt = optimoptions(optionsOpt,'PlotFcns', @optimplotfval);
optionsOpt = optimoptions(optionsOpt,'GradObj', 'off');
optionsOpt = optimoptions(optionsOpt,'Hessian', 'off');
optionsOpt = optimoptions(optionsOpt,'Algorithm', 'quasi-newton');

optionsOpt = optimoptions(optionsOpt, 'MaxFunEvals',300000);
optionsOpt = optimoptions(optionsOpt, 'TolFun', 1e-6);
optionsOpt = optimoptions(optionsOpt, 'TolX', 1e-6);
if ~isempty(maxIter) 
    optionsOpt = optimoptions(optionsOpt, 'maxIter', maxIter);
end

%% Optimization
F = @(coefs) rigidAlignmentDist( coefs(1:noCurves-1), ...
    coefs(noCurves:2*(noCurves-1)), ...
    reshape(coefs(2*(noCurves-1)+1:end), [noCurves-1, dSpace]), ...
    dAligned(1:noCurves), splineData, ...
    optTra, optRot, optShift, a );

init_coefs = zeros([(2+dSpace)*(noCurves-1), 1]);

% Deal with global rotations
if optRot && globalRot
    gaInit = findInitRot( dList, splineData, a, options);
    if optShift
        init_coefs(1) = gaInit.alpha;
    end
    init_coefs(noCurves-1+1) = gaInit.beta;
end;

optimal_coefs = fminunc( F, init_coefs, optionsOpt );

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
                                            splineData );
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
                                 dList, splineData, ...
                                 optTra, optRot, optShift, a )
% Only for dSpace=2 and periodic curves

dSpace = size(dList{1}, 2);
N = splineData.N;
noCurves = length(dList);

%% No transformation is applied to first curve
alpha = [0; alpha];
beta = [0; beta];
lambda = [ zeros([1, dSpace]); lambda];
                                 
%% Apply translations (lambda) and rotations (beta)
dTransformed = dList;
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

for jj = noCurves:-1:1
    if optShift
        dTransformed{jj} = ...
            curveApplyShift( dTransformed{jj}, alpha(jj), splineData );
    end
end

D = 0;
for jj = 1:noCurves-1
    for kk = jj+1:noCurves
        [E, ~] = curveFlatH2Norm( dTransformed{jj} - dTransformed{kk}, ...
                                  splineData, 'a', a );
        D = D + E;
    end
end

end

function gaInit = findInitRot( dList, splineData, a, options )
    d0 = dList{1};
    d1 = dList{2};
    
    options2 = options;
    options2.rigidGlobalRot = false;
    options2.maxIter = 40;
    
    noRotGuess = 6;
    
    distList = zeros([noRotGuess, 1]);
    for jj = 1:noRotGuess
        angle = 2*pi / noRotGuess * (jj-1);
        
        rotation = [  cos(angle), sin(angle); ...
                     -sin(angle), cos(angle) ];
        dTmp = d1* rotation;
        
        if isfield(options, 'optShift') && options.optShift
            dTmp = curveApplyShift(dTmp, -angle, splineData);
        end
        
        [dTmpList, ~] = rigidAlignment({d0, dTmp}, splineData, ...
            'options', options2);
        
        distList(jj) = curveFlatH2Norm( dTmpList{1} - dTmpList{2}, ...
                                        splineData, 'a', a );
    end
    
    [~, ind] = min(distList);
    
    gaInit = struct( 'phi', [], 'beta', [], 'v', [], 'alpha', []);
    gaInit.beta = 2*pi / noRotGuess * (ind-1);
    if isfield(options, 'optShift') && options.optShift
        gaInit.alpha = -gaInit.beta;
    end
end

function [dAligned, gaOpt] = shiftOnlyTwoCurves(dList, splineData, a)
    d0 = dList{1};
    d1 = dList{2};
    
    for jj = splineData.N:-1:1
        dist(jj) = curveFlatH2Norm(d0 - circshift(d1,jj-1,1), ...
                        splineData, 'a', a);
    end
    [~, optShift] = min(dist);
    alpha = optShift * 2*pi / splineData.N;
    
    dAligned = {d0, circshift(d1, optShift-1, 1)};
    gaOpt = { struct( 'phi', [], 'beta', [], 'v', [], 'alpha', []), ...
              struct( 'phi', [], 'beta', [], 'v', [], 'alpha', alpha) };
    
end
