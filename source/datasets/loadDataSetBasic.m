%% loadDataSet
% Loads some basic examples.
%
% Optional arguments
%   'curves' = {'H', 'U'} or an array containing a list of strings.
%
function dList = loadDataSetBasic( splineData, varargin )

%% Optional arguments
constSpeed = false; % Reparametrize to constant speed
curveSet = {};
returnCell = true;

allCurves = {'H', 'O', 'T', 'U', 'prop0', 'prop3', 'prop4', ...
             'circle', 'wrap', 'wrap_open'};
noise = 0;

ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'constspeed'
                constSpeed = true;
            case 'curves'
                ii = ii + 1;
                curveSet = varargin{ii};
            case 'noise'
                ii = ii + 1;
                noise = varargin{ii};
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    end
    ii = ii + 1;
end

if isempty(curveSet)
    curveSet = allCurves;
elseif ~iscell(curveSet)
    curveSet = {curveSet};
    returnCell = false;
end

dList = {};
for kk=length(curveSet):-1:1
    switch curveSet{kk}
        case {'H', 'O', 'T', 'U'}
            d0 = loadCurveLetters(curveSet{kk}, splineData, constSpeed);
        case {'prop0', 'prop3', 'prop4'}
            d0 = loadCurvePropeller(curveSet{kk}, noise, true, ...
                                    splineData, constSpeed);
        case {'circle', 'wrap'}
            d0 = loadCurveWrap(curveSet{kk}, splineData, constSpeed);
        case {'wrap_open'}
            d0 = loadCurveWrapOpen(splineData, constSpeed);
    end
    
    dList{kk} = d0;
end

if ~returnCell
    dList = dList{1};
end

end

function d = loadCurveLetters(curveCode, splineData, constSpeed)

noSub = 20;
f = [];

switch(curveCode)
    case 'H'
        t = linspace(-pi/2, 0, 3*noSub+1)'; t = t(1:end-1);
        f = [2*cos(t), 3 + 2*sin(t)];
        t = linspace(3, 7, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [2+0*t, t]];
        t = linspace(pi, 0, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [4+2*cos(t), 7+2*sin(t)]];
        t = linspace(7, 3, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [6+0*t, t]];
        t = linspace(3, -3, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [6+0*t, t]];
        t = linspace(-3, -7, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [6+0*t, t]];
        t = linspace(0, -pi, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [4+2*cos(t), -7+2*sin(t)]];
        t = linspace(-7, -3, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [2+0*t, t]];
        t = linspace(0, pi, 2*noSub+1)'; t = t(1:end-1);
        f = [f; 2*cos(t), -3 + 2*sin(t)];
        t = linspace(-3, -7, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [-2+0*t, t]];
        t = linspace(0, -pi, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [-4+2*cos(t), -7+2*sin(t)]];
        t = linspace(-7, -3, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [-6+0*t, t]];
        t = linspace(-3, 3, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [-6+0*t, t]];
        t = linspace(3, 7, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [-6+0*t, t]];
        t = linspace(pi, 0, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [-4+2*cos(t), 7+2*sin(t)]];
        t = linspace(7, 3, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [-2+0*t, t]];
        t = linspace(-pi, -pi/2, 3*noSub+1)'; t = t(1:end-1);
        f = [f; 2*cos(t), 3 + 2*sin(t)];
    case 'O'
        t = linspace(pi/2, -3*pi/2, 36*noSub+1)'; t = t(1:end-1);
        f = [6*cos(t), 9*sin(t)];
    case 'T'
        t = linspace(0, 4, 5*noSub+1)'; t=t(1:end-1);
        f = [t, 9+0*t];
        t = linspace(pi/2, -pi/2, 3*noSub+1)'; t = t(1:end-1);
        f = [f; [4+2*cos(t), 7+2*sin(t)]];
        t = linspace(pi/2, pi, noSub+1)'; t = t(1:end-1);
        f = [f; [4+2*cos(t), 3+2*sin(t)]];
        t = linspace(3, -3, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [2+0*t, t]];
        t = linspace(-3, -7, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [2+0*t, t]];
        t = linspace(0, -pi, 10*noSub+1)'; t = t(1:end-1);
        f = [f; [2*cos(t), -7+2*sin(t)]];
        t = linspace(-7, -3, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [-2+0*t, t]];
        t = linspace(-3, 3, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [-2+0*t, t]];
        t = linspace(0, pi/2, noSub+1)'; t = t(1:end-1);
        f = [f; [-4+2*cos(t), 3+2*sin(t)]];
        t = linspace(-pi/2, -3*pi/2, 3*noSub+1)'; t = t(1:end-1);
        f = [f; [-4+2*cos(t), 7+2*sin(t)]];
        t = linspace(-4, 0, 5*noSub+1)'; t=t(1:end-1);
        f = [f; [t, 9+0*t]];
    case 'U'
        t = linspace(-pi/2, 0, noSub+1)'; t = t(1:end-1);
        f = [2*cos(t), -3 + 2*sin(t)];
        t = linspace(-3, 7, 4*noSub+1)'; t=t(1:end-1);
        f = [f; [2+0*t, t]];
        t = linspace(pi, 0, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [4+2*cos(t), 7+2*sin(t)]];
        t = linspace(7, 3, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [6+0*t, t]];
        t = linspace(3, -3, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [6+0*t, t]];
        t = linspace(0, -pi, 14*noSub+1)'; t=t(1:end-1);
        f = [f; [6*cos(t), -3+6*sin(t)]];
        t = linspace(-3, 3, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [-6+0*t, t]];
        t = linspace(3, 7, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [-6+0*t, t]];
        t = linspace(pi, 0, 2*noSub+1)'; t=t(1:end-1);
        f = [f; [-4+2*cos(t), 7+2*sin(t)]];
        t = linspace(7, -3, 4*noSub+1)'; t=t(1:end-1);
        f = [f; [-2+0*t, t]];
        t = linspace(-pi, -pi/2, noSub+1)'; t = t(1:end-1);
        f = [f; 2*cos(t), -3 + 2*sin(t)];
end

a = 9;
b = 6;
ecc = sqrt(1 - (b/a)^2);
[~, len] = ellipke(ecc^2); 
len = 4*a*len; % Length of letter 'O'

d = constructSplineApproximation(f, splineData);
d = d / len * 2*pi;

if constSpeed
    d = curveReparamConstSpeed(d, splineData);
end

end

function d = loadCurvePropeller(curveCode, noise, rescale, ...
                                splineData, constSpeed)

t = linspace(0, 2*pi, 400)';
                            
switch(curveCode)
    case 'prop0'
        c = [ cos(t) .* (1 + noise * sin(9*t)), ...
              sin(t) .* (1 + noise * sin(9*t)) ];
          
    case 'prop3'
        c = [ cos(t) .* (1 - 0.3*sin(3*t)), ...
              sin(t) .* (1 - 0.3*sin(3*t)) ];
        n = [ -cos(t).*(1-0.3*sin(3*t)) - sin(t).*(-0.3*cos(3*t)*3), ...
              -sin(t).*(1-0.3*sin(3*t)) + cos(t).*(-0.3*cos(3*t)*3) ];
        n = n ./ (sqrt(sum(n.^2, 2)) * [1 1]);
        c = c - n .* ( noise*(sin(9*t))*[1 1] );
    case 'prop4'
        c = [ cos(t) .* (1 + 0.2*cos(4*t)), ...
              sin(t) .* (1 + 0.2*cos(4*t)) ];
        n = [ -cos(t).*(1+0.2*cos(4*t)) - sin(t).*(-0.2*sin(4*t)*4), ...
              -sin(t).*(1+0.2*cos(4*t)) + cos(t).*(-0.2*sin(4*t)*4) ];
        n = n ./ (sqrt(sum(n.^2, 2)) * [1 1]);
        c = c - n .* ( noise*(cos(12*t))*[1 1] );
end

d = constructSplineApproximation(c, splineData);

if rescale && abs(noise > eps)
    dNoNoise = loadCurvePropeller( curveCode, 0, false, ...
                                   splineData, constSpeed );
    dLength = curveLength(dNoNoise, splineData);
elseif rescale
    dLength = curveLength(d, splineData);
else
    dLength = 2*pi;
end

d = d / dLength * 2*pi;

if constSpeed
    d = curveReparamConstSpeed(d, splineData);
end

end

function d = loadCurveWrap(curveCode, splineData, constSpeed)

switch(curveCode)
    case 'wrap'
        f = @(t) 1/100*[ 17 - 34*cos(t) - 79*cos(2*t) + ...
            7*cos(3*t) - 2*cos(4*t) + 3*cos(5*t) - 24*sin(t) ...
            - 86*sin(2*t) + 13*sin(3*t) - sin(4*t) - 5*sin(5*t), ...
            -5 - 58*cos(t) - 53*cos(2*t) - 3*cos(3*t) + ...
            13*cos(4*t) - 11*cos(5*t) + 79*sin(t) + 8*sin(2*t) + ...
            19*sin(3*t) - 8*sin(4*t) ];

    case 'circle'
        f = @(t) 1/100*[150*cos(-t - 3*pi/4), ...
                        150*sin(-t - 3*pi/4)] ;
end

d = constructSplineApproximation(f, splineData);
lenCircle = 3*pi;

d = d / lenCircle * 2*pi;

if constSpeed
    d = curveReparamConstSpeed(d, splineData);
end

end

function d = loadCurveWrapOpen(splineData, constSpeed)

coefs = [-5.9927, -1.5533; -5.8487, -1.3981; -5.5024, -1.0943; ...
         -5.0188, -0.7303; -4.5122, -0.3979; -4.0164, -0.1287; ...
         -3.5013,  0.1130; -3.0229,  0.3091; -2.4752,  0.4510; ...
         -2.0589,  0.5041; -1.3852,  0.5365; -1.2942, -0.1100; ...
         -2.0106, -0.3150; -2.4127, -0.5394; -2.9607, -0.7326; ...
         -3.4515, -0.9488; -3.9364, -1.2342; -4.4921, -1.4958; ...
         -4.8903, -2.1448; -4.1443, -2.1877; -3.7509, -2.0268;...
         -3.1930, -1.8175; -2.7194, -1.5820; -2.2027, -1.3426;...
         -1.7060, -1.0874; -1.2035, -0.8389; -0.7143, -0.5805;...
         -0.2030, -0.3509;  0.1257, -0.2342;  0.2964, -0.2089 ];
sdTmp = constructEmptySplineData;
sdTmp.curveClosed = 0;
sdTmp.N = 30;
sdTmp.nS = 3;
sdTmp.quadDegree = [6 4];
sdTmp = constructKnots(sdTmp);
sdTmp = setupQuadData(sdTmp);

d = curveSpline2Spline(coefs, sdTmp, splineData);

if constSpeed
    d = curveReparamConstSpeed(d, splineData);
end

end