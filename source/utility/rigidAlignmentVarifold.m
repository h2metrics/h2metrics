%% rigidAlignmentVarifold
%
% Computes the rigidAlignment between a set of curves using the
% VarifoldDistance
%
% NOTE: Code works only for dSpace=2. In more than two dimensions rotations have
% to be rethought.
%
% Input
%   dList
%       List of curves
%   splineData
%       General information about the splines used.
%
% Output
%   dAligned
%       List of aligned curves
%   gaOpt
%       List of gamma structures to attain alignment
% Standard presets for optimization: Shifts, Rotations and Translations are
% minimized.
% Note: uses L^2 for shifts.
 
function [dAligned, gaOpt] = rigidAlignmentVarifold( dList, splineData, varargin) 

a= [1 0 0];
optTra = true;
optRot = true;
optShift= false;
maxIter = [];
globalRot = false;

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
if isfield(options, 'rigidMaxIter')
    maxIter = options.rigidMaxIter;
end

if isfield(options, 'rigidGlobalRot')
    globalRot = options.rigidGlobalRot;
end

dSpace = splineData.dSpace;
display='off';
N = splineData.N;

noCurves = length(dList);
[~,center] = curveCenter( dList{1}, splineData );
for jj = 2:noCurves
   [dList{jj},~] = curveCenter( dList{jj}, splineData );
   dList{jj}(:,1) = dList{jj}(:,1)+center(1);
   dList{jj}(:,2) = dList{jj}(:,2)+center(2);
end
dAligned = dList;
gaOpt = cell([noCurves, 1]);

if noCurves < 2
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
optionsOpt = optimoptions(optionsOpt, 'TolFun', 1e-12);
optionsOpt = optimoptions(optionsOpt, 'TolX', 1e-12);
if ~isempty(maxIter) 
    optionsOpt = optimoptions(optionsOpt, 'maxIter', maxIter);
end





%% Optimization
F = @(coefs) rigidAlignmentDist(coefs(noCurves:2*(noCurves-1)), ...
    reshape(coefs(2*(noCurves-1)+1:end), [noCurves-1, dSpace]), ...
    dAligned(1:noCurves), splineData, ...
    optTra, optRot);

init_coefs = zeros([(2+dSpace)*(noCurves-1), 1]);

% Deal with global rotations
if optRot && globalRot
    gaInit = findInitRot( dList, splineData, options);
    init_coefs(noCurves-1+1) = gaInit.beta;
end;



optimal_coefs = fminunc( F, init_coefs, optionsOpt );

%% Apply transformations to curve
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
end

%% Deal with shifts using the L2 metric
for jj = noCurves:-1:1
    if optShift
       [dShifted, shift]= shiftOnlyTwoCurves({dAligned{1},dAligned{jj}}, splineData, a);
       gaOpt{jj}.alpha = shift; 
       dAligned{jj}= dShifted;
    end
end




end


function D = rigidAlignmentDist(beta, lambda, ...
                                 dList, splineData, ...
                                 optTra, optRot)
% Only for dSpace=2 and periodic curves

dSpace = size(dList{1}, 2);
N = splineData.N;
noCurves = length(dList);

%% No transformation is applied to first curve
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


D = 0;
for jj = 1:noCurves-1
    for kk = jj+1:noCurves
[E]  = varifoldDistanceSquared(dTransformed{jj}, dTransformed{kk}, splineData);
D = D + E;
    end
end

end

function gaInit = findInitRot( dList, splineData, options )
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
        
        distList(jj) = varifoldDistanceSquared(dTmpList{1}, dTmpList{2}, splineData);
    end
    
    [~, ind] = min(distList);
    
    gaInit = struct( 'phi', [], 'beta', [], 'v', [], 'alpha', []);
    gaInit.beta = 2*pi / noRotGuess * (ind-1);
    if isfield(options, 'optShift') && options.optShift
        gaInit.alpha = -gaInit.beta;
    end
end

function [dShifted, alpha] = shiftOnlyTwoCurves(dList, splineData, a)
    d0 = dList{1};
    d1 = dList{2};
    
    for jj = splineData.N:-1:1
        dist(jj) = curveFlatH2Norm(d0 - circshift(d1,jj-1,1), ...
                        splineData, 'a', a);
    end
    [~, optShift] = min(dist);
    alpha = (optShift-1) * 2*pi / splineData.N;
    
    dShifted = circshift(d1, optShift-1, 1);    
end



