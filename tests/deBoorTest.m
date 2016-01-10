%% Test 1 - periodic splines
N = 20;
nS = 5;
knotsS =  [(-nS):(N+nS)]/N*2*pi;

d_per = rand([N,1]);
d_nonper = [ d_per; d_per(1:nS) ];

noPts = 1000;
t = linspace(0, 2*pi, noPts);

y1 = deBoor(knotsS, nS, d_per, t, 4, 'periodic', true);
B = spcol(knotsS, nS+1, brk2knt( t, 4 ));
y2 = B * d_nonper;

assert( norm(y1-y2) < 1e-14 );

%% Test 2 - nonperiodic
Nt = 100; %Number of time control points
nT = 5; %time degree
knotsT = [ zeros(1, nT), linspace(0,1,Nt - nT + 1), ones(1, nT)];

noPts = 100;
t = linspace(0, 1, noPts);
d_time = rand([Nt,1]);

y1 = deBoor(knotsT, nT, d_time, t, 4, 'periodic', false);
B = spcol(knotsT, nT+1, brk2knt( t, 4 ));
y2 = B * d_time;

assert( norm(y1-y2) < 1e-14 );

%% Test 3 - spcol with nonincreasing points

Nt = 100; %Number of time control points
nT = 5; %time degree
knotsT = [ zeros(1, nT), linspace(0,1,Nt - nT + 1), ones(1, nT)];

noPts = 100;
t = rand([1, noPts]);
d_time = rand([Nt,2]);
order = 2;

[t_sort, I] = sort(t);

y1_sort = deBoor(knotsT, nT, d_time, t, order, ...
                 'periodic', false, 'method', 'spcol');
y1 = zeros([noPts*order, 2]);
for jj = 1:order
    y1(order*(I-1)+jj,:) = y1_sort(jj:order:end,:);
end

B = spcol(knotsT, nT+1, brk2knt( t_sort, order ));
y2 = B * d_time;
        
assert( norm(y1-y2) < 1e-14 );