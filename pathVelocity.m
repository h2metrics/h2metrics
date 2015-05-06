function [ d ] = pathVelocity(dPath,evalT, splineData )
%Extract control points for the velocity field of the path at each time t.

%TODO: determine if dPath contains a diffeomorphism

%d = zeros( size(dPath,1),size(dPath,2), length(t));

controlPointWeights = spcol( splineData.knotsT, splineData.nT+1, evalT);

d_x = reshape( dPath(:,1), splineData.N,splineData.Nt)*controlPointWeights';
d_y = reshape( dPath(:,2), splineData.N,splineData.Nt)*controlPointWeights';

d(:,:,1) = d_x;
d(:,:,2) = d_y;

d = permute(d, [1 3 2]);

end

