function [innerKnotsS, innerKnotsT, quadData] = setupQuadData( spData )
%TODO: Summary
%For optimization
% Input: d0,d1,m,mT,n,nT,quadDegree
% Output:
% innerKnotsS, innerKnotsT, quadData struct, (knotS, knotsT)

%Some code for handling optional inputs
% i = 1;
% while i <= length(varargin)
%     if (isa(varargin{i},'char') && isa(varargin{i+1},'double') && isa(varargin{i+2},'double'))
%          switch (lower(varargin{i}))
%              case 'eval'
%                  evalS = varargin{i+1};
%                  evalT = varargin{i+2};
%              otherwise
%                  error('Invalid option: ''%s''.',varargin{i});
%          end     
%     end
%     i = i + 1;  
% end

nS = spData.nS;
nT = spData.nT;
N = spData.N;
Nt_inner = spData.Nt_inner;
quadDegree = spData.quadDegree;

Nt = Nt_inner + 2;

%TODO: replace with spData
knotsS = spData.knotsS;
knotsT = spData.knotsT;
innerKnotsS = knotsS(nS+1:end-nS);
innerKnotsT = knotsT(nT+1:end-nT);

disp('setupQuadData.m, calling gaussianQuadratureData');
[quadPointsS, quadWeightsS] = gaussianQuadratureData( unique(innerKnotsS), 'degree', quadDegree(1));
[quadPointsT, quadWeightsT] = gaussianQuadratureData( unique(innerKnotsT), 'degree', quadDegree(2));
noQuadPointsS = length(quadPointsS);
noQuadPointsT = length(quadPointsT);

knots = {knotsS, knotsT};
%sp = spmak( knots, d_nonper) ; %Not used

%0,1,2 second order derivatives in space, 0-1 in time
B_S_quad = spcol( knotsS, nS+1, brk2knt( quadPointsS, 3 ),'sparse');
B_T_quad = spcol( knotsT, nT+1, brk2knt( quadPointsT, 2 ),'sparse');

%Periodic B-splines
B_S_quad_per = [B_S_quad(:,1:nS) + B_S_quad(:,end-nS+1:end), B_S_quad(:,nS+1:end-nS)];

%clear quadData;
quadData = struct('weights',[],'points',[],'B',[],'Bu',[],'Bt',[],'Buu',[],...
    'But',[],'Buut',[],'degree','knots','noControlPoints',[]);
%quadData = quadraturePointsAndWeights(annulusRefined, 'degree', [10,6]);
% plotData = struct('points',[],'B',[],'Bu',[],'Bt',[],'Buu',[],...
%     'But',[],'Buut',[]);

%quadWeightsBT = quadWeightsS'*quadWeightsT;

%quadData.weights = quadWeightsBT(:);
quadData.weights = {quadWeightsS,quadWeightsT};
quadData.points = {quadPointsS, quadPointsT};
%quadData.degree = quadDegree;
quadData.degree = [nS,nT];
quadData.knots = knots;
quadData.noControlPoints = [N,Nt];
%plotData.points = {linspace(0,1,100), linspace(0,1,50) };

disp('setupQuadData.m, starting the double loop, which is really slow!');
tic
for ii = N:-1:1
    for jj = Nt_inner+2:-1:1;
     
    quadData.B(:,ii+(jj-1)*N) = reshape(B_S_quad_per(1:3:end,ii)*B_T_quad(1:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    quadData.Bu(:,ii+(jj-1)*N) = reshape(B_S_quad_per(2:3:end,ii)*B_T_quad(1:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1)' ;
    quadData.Bt(:,ii+(jj-1)*N) = reshape(B_S_quad_per(1:3:end,ii)*B_T_quad(2:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1)' ;
    quadData.Buu(:,ii+(jj-1)*N) = reshape(B_S_quad_per(3:3:end,ii)*B_T_quad(1:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1)' ;
    quadData.But(:,ii+(jj-1)*N) = reshape(B_S_quad_per(2:3:end,ii)*B_T_quad(2:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1)' ;
    quadData.Buut(:,ii+(jj-1)*N) = reshape(B_S_quad_per(3:3:end,ii)*B_T_quad(2:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1)' ;
    end
end
toc

end

