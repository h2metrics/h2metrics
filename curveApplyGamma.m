function c = curveApplyGamma(d, gamma, splineData, quadData, varargin)

applyDiff = true;
applyShift = true;

% Some code for handling optional inputs
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'applydiff'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    applyDiff = logical(varargin{ii});
                else
                    error('Invalid value for option ''applyDiff''.');
                end
            case 'applyshift'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    applyShift = logical(varargin{ii});
                else
                    error('Invalid value for option ''applyShift''.');
                end
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1;  
    end
end

c = d;
if applyDiff && ~isempty(gamma.phi)
    if applyShift && ~isempy(gamma.alpha)
        c = curveComposeDiff( c, gamma.phi - gamma.alpha, ...
                              splineData, quadData );
    else
        c = curveComposeDiff( c, gamma.phi, ...
                              splineData, quadData );
    end
elseif applyShift && ~isempty(gamma.alpha)
    c = curveApplyShift(c, gamma.alpha, splineData, quadData);
end
if ~isempty(gamma.v)
    c = c + ones([splineData.N, 1]) * gamma.v';
end
if ~isempty(gamma.beta)
    rotation = [ cos(gamma.beta), sin(gamma.beta); ...
                 -sin(gamma.beta), cos(gamma.beta) ];
    c = c * rotation;
end
