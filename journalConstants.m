%% Section 1: Setup Circle to wrap
dataName = 'journalConstants'; %name for the result files
%% Load curves
load('kimiaShapesN60');
extractnames = fieldnames((kimiaControlPoints));
for ii = 1:length(extractnames)
    dShapes{ii} = kimiaControlPoints.(extractnames{ii});
end
noCurves = length(dShapes);
splineData = constructEmptySplineData;
splineData.N = 60; %no. control points, must be bigger than n+1
splineData.Nt = 20 + 2; %Number of time control points
splineData.Nphi = 20; %No. control points for diffeomorphisms
splineData.nS = 3; %spacial degree
splineData.nT = 2; %time degree
splineData.nPhi = 3; %diffemorphism degree
splineData.quadDegree = [8,4]; %Quadrature precission
splineData.dSpace = 2;
splineData.noInterpolS = 2 * splineData.N; % For composition
splineData = constructKnots(splineData);
dSpace = 2;
splineData = setupQuadData(splineData);
splineData.a = [1 1 1];
dShapes(30) = [];
extractnames(30)=[];
n = length(dShapes);
options = struct( 'optDiff', true, ...
                  'optTra', true, ...
                  'optRot', true, ...
                  'optShift', true, ...
                  'tolFun', 1e-12, ...
                  'tolX', 1e-12, ...
                  'display', 'iter-detailed', ... % 'iter-detailed'
                  'maxIter', 1000,'rigidA' ,[1, 0, 0]);
%% Unparam
d0=dShapes{17};
d1=dShapes{32};
A=1;
B=0;
C=1;
splineData.a = [A B C];

%%
splineData.a = [A B 10*C]; 
dinitpath = linearPathCircle(d0, d1, splineData);
[optE, optP1,optGa1] = geodesicBvpDiff(d0,d1,splineData,'options', options,'initpath',dinitpath );
plotPathSingleThickThin(optP1,splineData,5);
L_1 = sqrt(optE);
figname = [plotDir, 'FishToolA1B0C10.pdf'];
export_fig(figname) 

    


%%
splineData.a = [A B 100*C]; 
dinitpath = linearPathCircle(d0, d1, splineData);
[optE, optP1,optGa1] = geodesicBvpDiff(d0,d1,splineData,'options', options,'initpath',dinitpath );
plotPathSingleThickThin(optP1,splineData,5);
L_2 = sqrt(optE);
figname = [plotDir, 'FishToolA1B0C100.pdf'];
export_fig(figname) 
    

%%
splineData.a = [A B 1000*C]; 
dinitpath = linearPathCircle(d0, d1, splineData);
[optE, optP1,optGa1] = geodesicBvpDiff(d0,d1,splineData,'options', options,'initpath',dinitpath );
plotPathSingleThickThin(optP1,splineData,5);
L_3 = sqrt(optE);
figname = [plotDir, 'FishToolA1B0C1000.pdf'];
export_fig(figname)     


%%
splineData.a = [A B 10000*C]; 
dinitpath = linearPathCircle(d0, d1, splineData);
[optE, optP1,optGa1] = geodesicBvpDiff(d0,d1,splineData,'options', options,'initpath',dinitpath );
plotPathSingleThickThin(optP1,splineData,5);
L_4 = sqrt(optE);
figname = [plotDir, 'FishToolA1B0C10000.pdf'];
export_fig(figname)

    




