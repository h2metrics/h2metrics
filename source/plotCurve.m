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

if nargin < 3
    lS = 'k-';
end

if isa(d,'cell')
    noCurves = length(d);
    for ii=1:noCurves;
        nS = splineData.nS;
        knotsS = splineData.knotsS;
        lS = lineStyle{ii};
        % Plot parameters
        noPlotPtsS = 500;
        plotPtsS = linspace(0, 2*pi, noPlotPtsS+1);
        
        c0 = deBoor(knotsS, nS, d{ii}, plotPtsS, 1, 'periodic', true);
        pt0 = deBoor(knotsS, nS, d{ii}, 0, 1, 'periodic', true);
        
        %% Setup plotting
        hold on;
        axis equal;
        
        %% Do plotting
        plot(c0(:, 1), c0(:, 2), lS, 'LineWidth', 1);
        plot(pt0(1), pt0(2), 'ko', 'LineWidth', 1);
        
        hold off; 
    end    
else    

    nS = splineData.nS;
    knotsS = splineData.knotsS;
    
    % Plot parameters
    noPlotPtsS = 500;
    plotPtsS = linspace(0, 2*pi, noPlotPtsS+1);
    
    c0 = deBoor(knotsS, nS, d, plotPtsS, 1, 'periodic', true);
    pt0 = deBoor(knotsS, nS, d, 0, 1, 'periodic', true);
    
    %% Setup plotting
    hold on;
    axis equal;
    
    %% Do plotting
    plot(c0(:, 1), c0(:, 2), lS, 'LineWidth', 1);
    plot(pt0(1), pt0(2), 'ko', 'LineWidth', 1);
    
    hold off;
end   
