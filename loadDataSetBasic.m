%% loadDataSet
% Loads some basic examples.
%
% Optional arguments
%   'curves' = {'H', 'U'} or an array containing a list of strings.
%
function dList = loadDataSetBasic( splineData, ~, varargin )

%% Optional arguments
constSpeed = false; % Reparametrize to constant speed
curveSet = {};
returnCell = true;

allCurves = {'H', 'O', 'T', 'U', 'prop0', 'prop3', 'prop4'};
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
    d = curveReparamConstSpeed(d, splineData, splineData.quadData);
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
    dLength = curveLength(dNoNoise, splineData, splineData.quadData);
elseif rescale
    dLength = curveLength(d, splineData, splineData.quadData);
else
    dLength = 2*pi;
end

d = d / dLength * 2*pi;

if constSpeed
    d = curveReparamConstSpeed(d, splineData, splineData.quadData);
end

end
