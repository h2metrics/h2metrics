%% Load curves
load('kimiasplinesN40ns3');
extractnames = fieldnames(rmfield(kimiaControlPoints,{'N','nS'}));
for ii = 1:length(extractnames)
    dShapes{ii} = kimiaControlPoints.(extractnames{ii});
end
noCurves = length(dShapes);
splineData = constructEmptySplineData;
splineData.N = 40; %no. control points, must be bigger than n+1
splineData.Nt = 10 + 2; %Number of time control points
splineData.Nphi = 6; %No. control points for diffeomorphisms
splineData.nS = 3; %spacial degree
splineData.nT = 2; %time degree
splineData.nPhi = 3; %diffemorphism degree
splineData.quadDegree = [8,4]; %Quadrature precission
splineData.dSpace = 2;
splineData.noInterpolS = 2 * splineData.N; % For composition
splineData = constructKnots(splineData);
dSpace = 2;
[quadData, quadDataTensor] = setupQuadData(splineData);
splineData.a = [1 1 1];

%%
%CAT to COW (load curves and rotate)
% dCurves = {kimiaControlPoints.cat1,kimiaControlPoints.cow1};
% %Rotate and scale curves
% A=zeros(2);
% A(1,2)=1;
% A(2,1)=-1;
% noCurves = length(dCurves);
% for kk = noCurves:-1:1
%     cQuad = quadData.B_S * dCurves{kk};
%     cQuad_u = quadData.Bu_S * dCurves{kk};
%     cSpeed = sum( cQuad_u.^2 , 2).^(1/2); 
%     cLength = sum(cSpeed .* quadData.quadWeightsS);
%     dCurves{kk} = A*(dCurves{kk}');
%     dCurves{kk} = dCurves{kk}';
%     dCurves{kk} = dCurves{kk} / cLength;
% end
%Determine constants
% [A B C] = determineConstants(dCurves,splineData,quadData,quadDataTensor);
% splineData.a = [A B C];
% [A B C]
% [E_geo, dPath] = geodesicBvpAmpl(dCurves{1},dCurves{2},splineData,quadData,quadDataTensor);
% plotParticlePath(dPath,splineData,6,40);
% export_fig CatCowA1B1C1Particle.pdf
% %plotPathSingle(dPath,splineData,6)
% plotPathRow(dPath,splineData,6);
% export_fig CatCowA1B1C1Row.pdf


% 
% 
% %%
% %Fish to Bunny (load curves)
dCurves = {kimiaControlPoints.bonefishes,kimiaControlPoints.bunny04};
%scale curves
noCurves = length(dCurves);
for kk = noCurves:-1:1
    cQuad = quadData.B_S * dCurves{kk};
    cQuad_u = quadData.Bu_S * dCurves{kk};
    cSpeed = sum( cQuad_u.^2 , 2).^(1/2); 
    cLength = sum(cSpeed .* quadData.quadWeightsS);
    dCurves{kk} = dCurves{kk} / cLength;
end
%Determine constants
[A B C] = determineConstants(dCurves,splineData,quadData,quadDataTensor);
[A B C]
splineData.a = [1*A 1*B 99*C];
[E_geo, dPath] = geodesicBvpAmpl(dCurves{1},dCurves{2},splineData,quadData,quadDataTensor);
plotParticlePath(dPath,splineData,6,40);
export_fig FishBunnyA1C1C99Particle.pdf
plotPathRow(dPath,splineData,6);
export_fig FishBunnyA1B1C99Row.pdf












