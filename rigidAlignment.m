function dAligned = rigidAlignment(dList, splineData, quadData, varargin) 

useComp = false;
maxIter = [];
display = 'iter';

ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'usecomp'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    useComp = logical(varargin{ii});
                else
                    error('Invalid value for option ''useComp''.');
                end
            case 'maxiter'
                ii = ii + 1;
                if isnumeric(varargin{ii})
                    maxIter = varargin{ii};
                else
                    error('Invalid value for option ''maxIter''.');
                end
            case 'display'
                ii = ii + 1;
                display = varargin{ii};
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1; 
    end
end

dSpace = splineData.dSpace;
N = splineData.N;

noCurves = length(dList);
dAligned = dList;

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
    dAligned(1:noCurves), splineData, quadData);

init_coefs = zeros([(2+dSpace)*(noCurves-1), 1]);

optimal_coefs = fminunc( F, init_coefs, options );

%% Apply transformations to curve
alpha = [0; optimal_coefs(1:noCurves-1)];
beta = [0; optimal_coefs(noCurves:2*(noCurves-1))];
lambda = [0,0; reshape(optimal_coefs(2*(noCurves-1)+1:end), ...
               [noCurves-1, dSpace])];

for jj = 1:noCurves
    dAligned{jj} = dList{jj} + ...
        ones([N, 1]) * lambda(jj,:);
    rotation = [ cos(beta(jj)), -sin(beta(jj)); ...
                 sin(beta(jj)),  cos(beta(jj)) ];
    dAligned{jj} = dAligned{jj} * rotation;
end

if useComp
    % Shift by alpha via diffeomorphism
    splineData2 = splineData;
    splineData2.Nphi = 5;
    splineData2.nPhi = 3;
    splineData2.noInterpolS = 5 * max(splineData2.N, splineData2.Nphi);
    splineData2 = constructKnots(splineData2);
    quadData2 = setupQuadData(splineData2);

    for jj = 2:noCurves
        phi = ones([splineData2.Nphi, 1]) * alpha(jj);
        dAligned{jj} = curveComposeDiff( dAligned{jj}, phi, ...
                                         splineData2, quadData2);
    end
else
    % We assume uniform knots and simply shift the control point sequence
    shift = round(alpha * N / (2*pi));
    for jj = 1:noCurves
        dAligned{jj} = circshift(dAligned{jj}, ...
            [-shift(jj), 0]);
    end
end

end



function D = rigidAlignmentDist( alpha, beta, lambda, ...
                                     dList, splineData, quadData )
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
    dTransformed{jj} = dList{jj} + ones([N, 1]) * lambda(jj,:);
    rotation = [ cos(beta(jj)), -sin(beta(jj)); ...
                 sin(beta(jj)),  cos(beta(jj)) ];
    dTransformed{jj} = dTransformed{jj} * rotation;
end

%% Evaluate shifted (alpha) curves at quadrature sites
cList = {};
for jj = noCurves:-1:1
    cList{jj} = deBoor( knotsS, nS, dTransformed{jj}, ...
                        mod(quadPointsS + alpha(jj), 2*pi), 1, ...
                        'periodic', true );
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
