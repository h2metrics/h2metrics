%% Tests linearPath
%
% Output should be a diagonal

% inner knots are in interval [2, 3]
knotsT = [-1, -.2, 1.2, 2, 2.1, 2.2, 2.4, 2.8, 3, 3.2, 3.4, 3.4]';
nT = 3;
Nt = 8;

splineData = struct( 'nT', nT, 'Nt', Nt, 'knotsT', knotsT );

d0 = -1;
d1 = 2;

dPath = linearPath(d0, d1, splineData);

noPlotPts = 100;
t = linspace(2, 3, noPlotPts);
B = spcol(knotsT, nT+1, brk2knt(t, 1));
y = B * dPath;

figure(2)
plot(t, y);