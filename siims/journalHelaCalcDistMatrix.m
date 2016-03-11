function journalHelaCalcDistMatrix(results, bothDir)
                               
%%
indList = results.indList;
noCurves = results.noCurves;
splineData = results.splineData;
options = results.options;
dList = results.dList;

%% Calculate distance matrix
distMatrix = zeros(noCurves, noCurves);
timeMatrix = zeros(noCurves, noCurves);
infoMatrix = cell(noCurves, noCurves);
dPathMatrix = cell(noCurves, noCurves);
gammaMatrix = cell(noCurves, noCurves);

for jj = 1:noCurves
    for kk = 1:noCurves
        % Forget the diagonal
        if jj == kk
            continue;
        end
        
        if jj > kk && ~bothDir
            continue;
        end
        
        disp(['dist ', num2str(indList(jj)), ' to ', ...
              num2str(indList(kk))]);
        
        d1 = dList{jj};
        d2 = dList{kk};
        
        tic
        [optE, optPath, optGa, info] = geodesicBvp(d1, d2, ...
            splineData, 'options', options);
        time = toc;
        
        distMatrix(jj, kk) = sqrt(optE);
        timeMatrix(jj, kk) = time;
        infoMatrix{jj, kk} = info;
        dPathMatrix{jj, kk} = optPath;
        gammaMatrix{jj, kk} = optGa;

    end
end

results.distMatrix = distMatrix;
results.timeMatrix = timeMatrix;
results.infoMatrix = infoMatrix;
results.dPathMatrix = dPathMatrix;
results.gammaMatrix = gammaMatrix;
disp('Calculations finished');