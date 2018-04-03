%% linearPath
%
% Function computes the linear interpolation between d0 and d1 using the
% spline described by knotsT, Nt and nT.
%
% Input
%   d0, d1
%       Matrices to be linearly interpolated
%   splineData
%       Contains knotsT, Nt and nT
%
% Output
%   dPath
%       The interpolated path. The dimensions of dPath are
%           [size(d0,1)*Nt, size(d0,2)]

function dPath = linearPath( d0, d1, splineData )

if ~isequal( size(d0), size(d1) )
    error('Dimension mismatch.');
end

N = size(d0, 1);
Nt = splineData.Nt;
nT = splineData.nT;
knotsT = splineData.knotsT;

t_min = knotsT(nT+1); % Limits for inner knots
t_max = knotsT(end-nT);
t_diff = t_max - t_min;

d_greville = aveknt(knotsT, nT+1)'; % This can be replace by an explicit 
                                    % formula, if necessary for AMPL.

for jj = Nt:-1:1
    for kk = N:-1:1
        dPath((jj-1)*N+kk,:) = ...
            d0(kk,:) * (t_max - d_greville(jj)) / t_diff + ...
            d1(kk,:) * (d_greville(jj) - t_min) / t_diff;
    end
end

end