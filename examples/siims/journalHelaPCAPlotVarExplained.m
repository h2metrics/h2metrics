disp(mfilename);

results_1p = matfile([prefixDir, 'hela/hela_1p.mat']);
results_2p = matfile([prefixDir, 'hela/hela_2p.mat']);
results_1d = matfile([prefixDir, 'hela/hela_1d.mat']);
results_2d = matfile([prefixDir, 'hela/hela_2d.mat']);

%% Extract things
varExpl_1p = results_1p.meanVarExplained;
varExpl_2p = results_2p.meanVarExplained;
varExpl_1d = results_1d.meanVarExplained;
varExpl_2d = results_2d.meanVarExplained;

%% Plot 1, variance explained
plotDir = [prefixDir, 'final_plots/'];
noPC = 10;

% Setup plotting
lineWidth = 400;
figRelSize = 0.49;

figRatio = 4/3;
sx = figRelSize * lineWidth;
sy = sx / figRatio;
handle = figure( 'PaperUnits', 'points', 'PaperSize', [sx, sy], ...
                 'Units', 'points', 'Position', [0, 0, sx, sy], ...
                 'Color', 'white' );
handle.Visible = 'off';

hold on;
plot(varExpl_1d(1:noPC), 'k-x');
plot(varExpl_2d(1:noPC), 'k--x');
plot(varExpl_1p(1:noPC), 'k-o');
plot(varExpl_2p(1:noPC), 'k--o');


legend({'$a_2=2^{-12}$, mod Diff', ...
        '$a_2=2^{-8}$, mod Diff', ...
        '$a_2=2^{-12}$, param', ...
        '$a_2=2^{-8}$, param'}, ...
        'Location', 'best', 'Interpreter', 'latex');

% Finish the plot
hold off;

figname = [ plotDir, 'pca_variance_explained.eps'];
export_fig(figname);