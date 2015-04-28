function d = constructSplineApproximation(f,splineData,varargin);
% Given a function f, construct an approximation of f using the data in
% splineData.
%
% Input: f, function handle
%        splineData, splineData struct
%
% Output: d, control points of the computed approximation

if isa(f, 'function_handle') %f is a handle
    % Interpolate data from parametrizations with splines
    noInterpolPoints = max(2*splineData.N,200);
    interpolS = linspace( 0, 2*pi, noInterpolPoints+1)';
    interpolS = interpolS(1:end-1); %remove last point, correponds to first point
    B_interpol = spcol( splineData.knotsS, splineData.nS+1, brk2knt( interpolS, 1 ),'sparse');
    B_interpol_p = [B_interpol(:,1:splineData.nS) + B_interpol(:,end-splineData.nS+1:end), B_interpol(:,splineData.nS+1:end-splineData.nS)];
    
    %Evaluate parametrization
    data = f(interpolS);
    
    %Solve the linear interpolation problem
    d = B_interpol_p\data;
    
elseif isa(f, 'numeric') % f is a list of points
    noInterpolPoints = size(f,1);
    interpolS = linspace( 0, 2*pi, noInterpolPoints+1)';
    interpolS = interpolS(1:end-1); %remove last point, correponds to first point
    B_interpol = spcol( splineData.knotsS, splineData.nS+1, brk2knt( interpolS, 1 ),'sparse');
    B_interpol_p = [B_interpol(:,1:splineData.nS) + B_interpol(:,end-splineData.nS+1:end), B_interpol(:,splineData.nS+1:end-splineData.nS)];
    
    d = B_interpol_p\f;
else
    error('Unknown f');
end





end

