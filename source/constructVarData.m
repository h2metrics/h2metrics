%% constructVarData
%
% Function constructs an varData struct with preset parameters. Some default parameters
% are being set; others are set to empty arrays.
%
% If splineData is provided, then pts and the connectivity matrix are
% constructed.
%
% Default parameters
%   kernelGeom = 'gaussian_oriented'
%   kernelGrass = 'binet'
%
% Output
%   varData
%       The created struct.
%
function [ varData ] = constructVarData(splineData, varargin)

% Handle optional inputs
p = inputParser;
addParameter(p, 'kernelGeom', 'gaussian');
addParameter(p, 'kernelGrass', 'gaussian_oriented');
addParameter(p, 'kernelSizeGeom', 0.1);
addParameter(p, 'kernelSizeGrass', 0.3);
addParameter(p, 'noPts', splineData.N);
addParameter(p, 'pts', []);
addParameter(p, 'G', []);
parse(p, varargin{:});

varData = p.Results;

% If noPts is not set, try setting with splineData
if isempty(varData.noPts) && ~isempty(splineData.N)
    if isempty(splineData.nS)
        nS = 3;
    else
        nS = splineData.nS;
    end
    varData.noPts = nS * splineData.N;
end

% If nothing still, return
if isempty(varData.noPts)
    return
end

% Set evaluation points and connectivity matrix
if splineData.curveClosed
    varData.pts = linspace(0, 2*pi, varData.noPts+1);
    varData.pts = varData.pts(1:end-1); % Remove the last repeated point

    % Connectivity matrix; note, last connected to first
    varData.G = [(1:varData.noPts)', circshift((1:varData.noPts)', -1)];
else
    varData.pts = linspace(0, 2*pi, varData.noPts);

    % Connectivity matrix
    varData.G = [(1:varData.noPts-1)', (2:varData.noPts)'];
end

end

