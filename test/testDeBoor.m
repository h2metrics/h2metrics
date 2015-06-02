%% testDeBoor
%
% Script shows how to use deBoor.

%% Setup variables
N = 20; %number of periodic control points, must be bigger than n+1
nS = 5; %spacial degree

Nt = 20; %Number of time control points
nT = 5; %time degree

knotsS =  [(-nS):(N+nS)]/N*2*pi; %normalize, domain of definition is [0,2*pi]x[0,2*pi]
knotsT = [ zeros(1, nT), linspace(0,1,Nt - nT + 1), ones(1, nT)];
% innerKnotsS = knotsS(nS+1:end-nS);

%% Given control points
d_per = rand([N,1]);
d_nonper = [ d_per; d_per(1:nS) ];
d_time = rand([Nt, 1]);

%% Where to evaluate
t = 2*pi;

%% Evaluate de Boor's algorithm and compare with spcol
r = floor(t*N/(2*pi)+nS)+1;
if r == N + nS + 1
    r = N + nS;
end
s = knotsS(r-nS+1:r+nS);
d_work = d_nonper(r-nS:r)';

for kk = 1:nS
    al(kk:nS) = (t - s(kk:nS)) ./ (s(nS+1:2*nS-kk+1) - s(kk:nS));
    d_work(kk+1:nS+1) = (1-al(kk:nS)) .* d_work(kk:nS) + al(kk:nS) .* d_work(kk+1:nS+1);
end

B = spcol(knotsS, nS+1, brk2knt( t, 1 ));
disp(B * d_nonper - d_work(nS+1));

%% Test function, periodic
noPts = 1000;
t = linspace(0, 2*pi, noPts);
d_per = rand([N,1]);
d_nonper = [ d_per; d_per(1:nS,:) ];

tic
y = deBoor(knotsS, nS, d_per, t, 4, 'periodic', true);
toc
tic
y2 = zeros([4*noPts,1]);
B = spcol(knotsS, nS+1, brk2knt( t, 4 ));
y2(4:4:end,:) = B(4:4:end,:) * d_nonper;
y2(3:4:end,:) = B(3:4:end,:) * d_nonper;
y2(2:4:end,:) = B(2:4:end,:) * d_nonper;
y2(1:4:end,:) = B(1:4:end,:) * d_nonper;
toc
disp(max(abs(y(1:4:end,:)-y2(1:4:end,:))))
disp(max(abs(y(2:4:end,:)-y2(2:4:end,:))))
disp(max(abs(y(3:4:end,:)-y2(3:4:end,:))))
disp(max(abs(y(4:4:end,:)-y2(4:4:end,:))))

%% Test function, nonperiodic
noPts = 1000;
t = linspace(0, 1, noPts);
d_time = rand([Nt,1]);

tic
y = deBoor(knotsT, nT, d_time, t, 4, 'periodic', false);
toc
tic
y2 = zeros([4*noPts,1]);
B = spcol(knotsT, nT+1, brk2knt( t, 4 ));
y2(4:4:end,:) = B(4:4:end,:) * d_time;
y2(3:4:end,:) = B(3:4:end,:) * d_time;
y2(2:4:end,:) = B(2:4:end,:) * d_time;
y2(1:4:end,:) = B(1:4:end,:) * d_time;
toc
disp(max(abs(y(1:4:end,:)-y2(1:4:end,:))))
disp(max(abs(y(2:4:end,:)-y2(2:4:end,:))))
disp(max(abs(y(3:4:end,:)-y2(3:4:end,:))))
disp(max(abs(y(4:4:end,:)-y2(4:4:end,:))))