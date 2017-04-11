disp(mfilename);

%% Load files
results_1d = matfile([prefixDir, 'hela/hela_1d.mat']);
results_2d = matfile([prefixDir, 'hela/hela_2d.mat']);

resultsList = {results_1d, results_2d};

%% Print information
for jj = 1:length(resultsList)
    results = resultsList{jj};
    splineData = results.splineData;
    options = results.options;
    
    % Describe splineData
    disp(' ');
    disp('-----------------------------------------------');
    disp(['N=', num2str(splineData.N), ', nS=', num2str(splineData.nS)]);
    disp(['Nt=', num2str(splineData.Nt), ', nT=', num2str(splineData.nT)]);
    disp(['Nphi=', num2str(splineData.Nphi), ...
          ', nPhi=', num2str(splineData.nPhi)]);
    quadPoints = ceil((splineData.quadDegree-1)/2);
    disp(['quadDegree=(', num2str(splineData.quadDegree(1)), ', ',...
          num2str(splineData.quadDegree(2)), '), ',...
          'quadPoints=(', num2str(quadPoints(1)), ', ', ...
          num2str(quadPoints(2)), ')']);
    disp(['a(1)=', num2str(splineData.a(1)), ...
          ', a(2)=', num2str(splineData.a(2)), ...
          ', a(2)=', num2str(splineData.a(3))]);
    disp(['optDiff=', num2str(options.optDiff), ...
          ', optTra=', num2str(options.optTra), ...
          ', optRot=', num2str(options.optRot), ...
          ', optShift=', num2str(options.optShift) ]);
    disp(' ');
    
    dPathMatrix = results.dPathMatrix;
    noCurves = results.noCurves;
    
    rhoL2 = zeros([noCurves, noCurves]);
    rhoH2 = zeros([noCurves, noCurves]);
    
    for kk = 1:noCurves
        for ll = kk+1:noCurves
            [E, comp] = pathRiemH2Energy(dPathMatrix{kk,ll}, splineData);
            rhoL2(kk,ll) = comp(1) / E;
            rhoH2(kk,ll) = comp(3) / E;
        end
    end
    
    rhoL2vec = squareform(rhoL2');
    rhoH2vec = squareform(rhoH2');
    
    disp(['a(2)=', num2str(splineData.a(3))]);
    disp(['rhoH2mean=', num2str(mean(rhoH2vec))]);
    disp(['rhoH2std=', num2str(std(rhoH2vec))]);
    disp(' ');
    
    varExplained = results.meanVarExplained;
    disp(['Variance explained=', num2str(varExplained(1:6,1)')]);
    
    
end

