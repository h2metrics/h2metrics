clear;

%% Common variables
dataDir = '~/diss/openprojects/h2_numerics/data/';
prefixDir = '~/diss/openprojects/h2_numerics/journal2/';
plotDir = [prefixDir, 'final_plots/'];

%% Figure 1, images of cells
clearvars -except dataDir prefixDir plotDir
batch('journalHelaPrintImages');

%% Figure 2, geodesics N=12, N=40
clearvars -except dataDir prefixDir plotDir
batch('journalHelaSomeGeodesics');

%% Print results
clearvars -except dataDir prefixDir plotDir
batch('journalHelaPrintResults');

%% Figure, variance explained
clearvars -except dataDir prefixDir plotDir
batch('journalHelaPCAPlotVarExplained');

%% PCA figures
clearvars -except dataDir prefixDir plotDir
batch('journalHelaPCAPlot');