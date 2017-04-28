%% constructStandardSplineData
% Function constructs a standard splineData struct with N=100 and Nt=6
% Possible varargins: different valus for [N,Nt], i.e., 
% constructStandardSplineData([N,Nt]) or constructStandardSplineData(N)
function [ splineData ] = constructStandardSplineData(varargin)
if isempty(varargin) 
   N=100;
   Nt=6;
elseif length(varargin{1})==1    
   N=varargin{1}; 
   Nt=6;
else
   N=varargin{1}(1); 
   Nt=varargin{1}(2);  
end




splineData = constructEmptySplineData;
splineData.N = N;
splineData.nS = 3;
splineData.Nt = Nt;
splineData.nT = 2;
splineData.quadDegree = [6, 4];
splineData.Nphi = 15; % 
splineData.a = [1 1e-1 1e-5];
splineData.nPhi = 3;
splineData.lambda=100;
varData = constructEmptyVarData(splineData);
varData.kernelSize = 0.1;
splineData.varData = varData;


splineData = constructKnots(splineData);
splineData = setupQuadData(splineData);





options = struct( 'optDiff', true, ...
                  'optTra', true, ...
                  'optRot', true, ...
                  'optShift', true, ...
                  'tolFun', 1e-12, ...
                  'tolX', 1e-12, ...
                  'display', 'iter-detailed', ... % 'off', 'iter-detailed'
                  'maxIter', 300);

options.rigidA = [1 0 0];
              
rigidoptions =  struct('optTra', true, ...
                  'optRot', false, ...
                  'optShift', true, ...
                  'tolFun', 1e-3, ...
                  'tolX', 1e-3, ...
                  'maxIter', 15 );    

rigidoptions.rigidA = [1 0 0];
splineData.options=options;
splineData.rigidoptions=rigidoptions;