%% constructSplineApproximation
%
% Given a function or a set of points, construct a spline approximating f
% using the data in splineData.
%
% Input
%   f
%       Can be a function handle or a set of points.
%   splineData
%       Information about the spline to be constructed.
%
% Output
%   d
%       Control points of the resulting spline.
%
function d = constructSplineApproximation(f,splineData,varargin)

if isa(f, 'function_handle') % f is a function handle
    interpolS = splineData.interpolS;
    B_interpol = splineData.quadData.B_interpolS;
    
    data = f(interpolS);
    d = B_interpol \ data;
    
elseif isa(f, 'numeric') % f is a list of points
    nS = splineData.nS;
    noInterpolPoints = size(f, 1);
    
    if splineData.curveClosed
        interpolS = linspace(0, 2*pi, noInterpolPoints+1)';
        interpolS = interpolS(1:end-1); % Last point correponds to first
        
        B_interpol = spcol( splineData.knotsS, nS+1, ...
                            brk2knt( interpolS, 1 ), 'sparse');
        B_interpol = ...
            [ B_interpol(:,1:nS) + B_interpol(:,end-nS+1:end), ...
              B_interpol(:,nS+1:end-nS) ];
    else
        interpolS = linspace(0, 2*pi, noInterpolPoints)';
        
        B_interpol = spcol( splineData.knotsS, nS+1, ...
                            brk2knt( interpolS, 1 ), 'sparse');
    end
    
    % Solve the linear interpolation problem
    d = B_interpol \ f;
    
else
    error('Unknown f');
    
end

end

