function journalHelaCalcKarcherMean(results)

noCurves = results.noCurves;
splineData = results.splineData;
options = results.options;
dList = results.dList;

N = splineData.N;
dSpace = splineData.dSpace;

meanInitInd = results.meanInitInd;
meanInit = dList{meanInitInd};

%% Calculate Karcher mean
disp('Computing the Karcher mean');
[dMean, info] = karcherMeanManopt( dList, ...
    splineData, 'options', options, 'meanInit', meanInit );
disp('Karcher mean computed');

%% Save results
results.mean = dMean;
results.meanInfo = info;

%% Save all paths
dPathList = info(end).dPathList;
gaList = info(end).gaList;

vList = {};
for jj = noCurves:-1:1
    vList{jj} = pathVelocity(dPathList{jj}, 0, splineData);
end

results.meanGammaList = gaList;
results.meanPathList = dPathList;
results.meanVelList = vList;

%% Do some PCA
G = metricMatrixH2(dMean, splineData);
rootG = sqrtm(full(G));

Sigma = zeros([N*dSpace, N*dSpace]);
for jj = noCurves:-1:1
    v = reshape(vList{jj}, [N*dSpace, 1]);
    Sigma = Sigma + v * v';
end

Sigma = 1./(noCurves-1) * rootG * Sigma * rootG;

[U, Lambda] = eig(Sigma);
Lambda = real(diag(Lambda));

varExplained = cumsum(Lambda) / sum(Lambda);
results.meanVarExplained = varExplained;

%% Calculate principal components
noPC = size(U, 1);
meanPCList = {};
for jj = noPC:-1:1
    v0 = rootG \ U(:,jj);
    v0 = sqrt(Lambda(jj)) * v0;
    v0 = reshape(v0, [N, dSpace]);
    meanPCList{jj} = v0;
end
results.meanPCList = meanPCList;

%% Calculate 2d coordinates
V = U(:,1:2);
pts2d = zeros([2, noCurves]);
for jj = noCurves:-1:1
    v = reshape(vList{jj}, [N*dSpace, 1]);
    pts2d(:,jj) = V' * rootG * v;
    pts2d(:,jj) = pts2d(:,jj) ./ sqrt(Lambda(1:2));
end
results.meanPts2d = pts2d;