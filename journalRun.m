%% Setup
clear;

dataDir = '~/diss/openprojects/h2_numerics/data/';
prefixDir = '~/diss/openprojects/h2_numerics/journal/';

%% Create basic shapes
journalBasic;

%% Fig. 1 - Tool to Fish

%% Fig. 2 - Plotting the shapes

%% Fig. 3 - Convergence of energy

%% Fig. 4 - Convergence of BVP

%% Fig. 5 - Continuity of distance function
clearvars -except dataDir prefixDir

batch('journalNoise');   % Fig 5L
batch('journalNoise2');  % Fig 5R

%% Fig. 6 - Symmetry
clearvars -except dataDir prefixDir

batch('journalSymmetryDiff'); % Fig 6L, 6R

%% Fig. 7 - Geodesic forward/backward

%% Fig. 8 - Compatibiliy IVP/BVP
clearvars -except dataDir prefixDir

batch('journalForward2');        % Fig 8L
batch('journalInitalVelocity');  % Fig 8R

%% Fig. 9 - Karcher mean propellers