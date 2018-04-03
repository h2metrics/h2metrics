%% constructEmptyVarData
%
% Function constructs an empty varData struct. Some default parameters
% are being set; others are set to empty arrays.
%
% If splineData is provided, then pts and the connectivity matrix are
% constructed.
%
% Default parameters
%   kernelGeom = 'gaussian'
%   kernelGrass = 'binet'
%
% Output
%   varData
%       The created struct.
%
function [ varData ] = constructEmptyVarData(splineData)

varData = struct( ...
    'kernelGeom', 'gaussian', ...
    'kernelGrass', 'binet', ...
    'kernelSize', [], ...
    'noPts', [], ...
    'pts', [], ...
    'G', [] );

if nargin > 0 && ~isempty(splineData.N)
    varData.noPts = splineData.N;
    
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

end

