%% Setup
clear;

dataDir = '../data/';
prefixDir = 'journal/';
plotDir = [prefixDir, 'final_plots/'];

%% Create basic shapes
journalBasic;

%% Fig. 1 - Tool to Fish
clearvars -except dataDir prefixDir plotDir
batch('journalConstants');  

%% Fig. 2 - Plotting the shapes
clearvars -except dataDir prefixDir plotDir
batch('journalShapesPlot');

%% Fig. 3 - Convergence of energy
clearvars -except dataDir prefixDir plotDir
batch('journalEnergySplineDegreeIndependence');   % Fig 3L 3R

%% Fig. 4 - Convergence of BVP
clearvars -except dataDir prefixDir plotDir
batch('journalConvPropellerFinal');

%% Fig. 5 - Continuity of distance function
clearvars -except dataDir prefixDir plotDir
batch('journalNoise');   % Fig 5L

clearvars -except dataDir prefixDir plotDir
batch('journalNoise2');  % Fig 5R

%% Fig. 6 - Symmetry
clearvars -except dataDir prefixDir plotDir
batch('journalSymmetryDiff'); % Fig 6L, 6R

%% Fig. 7 - Geodesic forward/backward
clearvars -except dataDir prefixDir plotDir
batch('journalForwardBackward'); % Fig 7T, 7B

%% Fig. 8 - Compatibiliy IVP/BVP
clearvars -except dataDir prefixDir plotDir
batch('journalForward2');        % Fig 8L

clearvars -except dataDir prefixDir plotDir
batch('journalInitalVelocity');  % Fig 8R

%% Fig. 9 - Karcher mean propellers
clearvars -except dataDir prefixDir plotDir
batch('journalKarcherMean');        % Fig 9
