function splineData = finishSetup(splineData)

% Construct varData
if isempty(splineData.varData)
%     varData = constructEmptyVarData(splineData);
    varData = constructVarData(splineData);
else
    varData = splineData.varData;
end


% Standard Parameter for varData
splineData.varData = varData;

% Construct knots
splineData = constructKnots(splineData);

% Construct collocation matrices
splineData = setupQuadData(splineData);


%%optimization parameter (should be included in options)
% splineData.eps_match = 1e-2; %Soft constraint for similarity term
% splineData.tau_final = 1e-3;

splineData.options = struct( 'optDiff', true, ...
                  'optTra', true, ...
                  'optRot', true, ...
                  'optScal', false, ...
                  'useVarifold', true, ...
                  'hansoNvec', 500, ...
                  'hansoMaxIt', 1500, ...
                  'hansoNormTol', 1e-3, ...
                  'hansoPrtLevel', 1,...
                  'useAugmentedLagrangian',true,...
                  'eps_match',1e-2,...
                  'tau_final', 1e-3);

