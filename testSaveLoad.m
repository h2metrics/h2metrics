%% Setup Variables

dSpace = 2;
N = 12; %no. control points, must be bigger than n+1
Nt_inner = 6; % Number of time control points between d0 and d1.
Nt = Nt_inner + 2; %Number of time control points
Nphi = 10; %Number of diffeo control points
nS = 5; %spacial degree
nT = 2; %time degree
nPhi = 5; %diffeo degree
quadDegree = [6,4]; %Quadrature precission
phi_eps = 1e-12;

splineData = struct( 'dSpace',dSpace, 'nS',nS, 'nT',nT, 'nPhi',nPhi,...
    'N',N, 'Nt',Nt, 'Nt_inner',Nt_inner, 'Nphi',Nphi, ...
    'quadDegree',quadDegree, 'knotsS',[], 'knotsT',[], 'knotsPhi', [],...
    'innerKnotsS',[], 'innerKnotsT',[], 'innerKnotsPhi',[], ...
    'quadData',[], 'phi_eps', phi_eps );
splineData.knotsS =  [(-nS):(N+nS)]/N*2*pi; %normalize, domain of definition is [0,2*pi]x[0,2*pi]
splineData.knotsT = [ zeros(1, nT), linspace(0,1,Nt_inner+2 - nT + 1), ones(1, nT)];
splineData.knotsPhi = [(-nPhi):(Nphi+nPhi)]/Nphi*2*pi;
splineData.innerKnotsS = splineData.knotsS(nS+1:end-nS);
splineData.innerKnotsT = splineData.knotsT(nT+1:end-nT);
splineData.innerKnotsPhi = splineData.knotsPhi(nPhi+1:end-nPhi);
splineData.Nt_inner = Nt_inner;

%% Create object to save
d0 = rand([N, dSpace]);
dPath = rand([N*Nt, dSpace]);
phiPath = rand([Nphi, dSpace]);

name1 = 'testfile1';
name2 = 'testfile2';
name3 = 'testfile3';
workdir = 'testdir';
dirExists = exist(workdir, 'dir');

if ~dirExists
    mkdir(workdir);
end

%% Save and load

saveCurve(name1, d0, splineData, 'workdir', workdir);
d1 = loadCurve(name1, 'workdir', workdir, 'autocomplete', true);
disp(max(max(abs(d0-d1))));

savePath(name2, dPath, splineData, 'workdir', workdir);
dPath1 = loadPath(name2, 'workdir', workdir, 'autocomplete', true);
disp(max(max(abs(dPath - dPath1))));

savePathDiff(name3, dPath, phiPath, splineData, 'workdir', workdir);
[dPath1, phiPath1, splineData] = loadPathDiff(name3, 'workdir', workdir, ...
                                              'autocomplete', true);
disp(max(max(abs(dPath - dPath1))));
disp(max(max(abs(phiPath - phiPath1))));
