%% plotCurve
%
% Helper function to plot a curve
%
% Input
%   d
%       A curve or a list of curves
%   splineData
%       General information about the splines used.
%
% Optional inputs
%   lineStyle = 'k-' (default)
%       lineStyle parameter to be passed to plot.
%
function plotCurve(d, splineData, lineStyle)

% Plot parameters
noPlotPtsS = 500;
plotPtsS = linspace(0, 2*pi, noPlotPtsS+1);
nS = splineData.nS;
knotsS = splineData.knotsS;

if ~isa(d,'cell')
   d={d}; 
   if nargin > 2
   lineStyle={lineStyle};
   end   
end        
noCurves = length(d);
for ii=1:noCurves
        if nargin < 3
            lS = 'k-';
        else
            lS = lineStyle{ii}; 
        end
        c0 = deBoor(knotsS, nS, d{ii}, plotPtsS, 1, 'periodic', true);
        pt0 = deBoor(knotsS, nS, d{ii}, 0, 1, 'periodic', true);
        %% Setup plotting
        hold on;
        axis equal;        
        %% Do plotting
        plot(c0(:, 1), c0(:, 2), lS);
        plot(pt0(1), pt0(2), 'ko');      
        hold off; 
end
end


