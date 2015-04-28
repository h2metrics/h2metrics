function [quadData, quadDataTensor] = setupQuadData( splineData)
%TODO: Summary
%For optimization
% Input: spData
% Output:
% innerKnotsS, innerKnotsT, quadData struct, (knotS, knotsT)

%clear quadData;
quadData = struct('quadPointsS',[],'quadPointsT',[],...
    'noQuadPointsS',[],'noQuadPointsT',[],...,
    'quadWeightsS',[],'quadWeightsT',[],...
    'B_S',[],'Bu_S',[],'Buu_S',[],'Buuu_S',[],...
    'B_T',[],'Bt_T',[],...
    'B_phi',[],'Bu_phi',[],'Buu_phi',[],'Buuu_phi',[]);

nS = splineData.nS;
nT = splineData.nT;
nPhi = splineData.nPhi;
N = splineData.N;
Nt = splineData.Nt;
Nphi = splineData.Nphi;
quadDegree = splineData.quadDegree;

%TODO: replace with spData
knotsS = splineData.knotsS;
knotsT = splineData.knotsT;
innerKnotsS = splineData.innerKnotsS;
innerKnotsT = splineData.innerKnotsT;

[quadPointsS, quadWeightsS] = gaussianQuadratureData( unique(innerKnotsS), 'degree', quadDegree(1));
[quadPointsT, quadWeightsT] = gaussianQuadratureData( unique(innerKnotsT), 'degree', quadDegree(2));
noQuadPointsS = length(quadPointsS);
noQuadPointsT = length(quadPointsT);

quadData.quadPointsS = quadPointsS';
quadData.quadPointsT = quadPointsT';
quadData.noQuadPointsS = noQuadPointsS;
quadData.noQuadPointsT = noQuadPointsT;
quadData.quadWeightsS = quadWeightsS';
quadData.quadWeightsT = quadWeightsT';

%0,1,2 second order derivatives in space, 0-1 in time
noSder = 3;
noPhider = 3;
noTder = 1;
B_S_quad = spcol( knotsS, nS+1, brk2knt( quadPointsS, noSder+1 ),'sparse');
B_Phi_quad = spcol( knotsS, nPhi+1, brk2knt( quadPointsS, noPhider+1 ),'sparse');
B_T_quad = spcol( knotsT, nT+1, brk2knt( quadPointsT, noTder+1 ),'sparse');
%Periodic B-splines
B_S_quad_per = [B_S_quad(:,1:nS) + B_S_quad(:,end-nS+1:end), B_S_quad(:,nS+1:end-nS)];
B_Phi_quad_per = [B_Phi_quad(:,1:nPhi) + B_Phi_quad(:,end-nPhi+1:end), B_Phi_quad(:,nPhi+1:end-nPhi)];

quadData.B_S = B_S_quad_per(1:4:end,:);
quadData.Bu_S = B_S_quad_per(2:4:end,:);
quadData.Buu_S = B_S_quad_per(3:4:end,:);
quadData.Buuu_S = B_S_quad_per(4:4:end,:);

quadData.B_T = B_T_quad(1:2:end,:);
quadData.Bt_T = B_T_quad(2:2:end,:);

quadData.B_phi = B_Phi_quad_per(1:4:end,:);
quadData.Bu_phi = B_Phi_quad_per(2:4:end,:);
quadData.Buu_phi = B_Phi_quad_per(3:4:end,:);
quadData.Buuu_phi = B_Phi_quad_per(4:4:end,:);

%quadData = quadraturePointsAndWeights(annulusRefined, 'degree', [10,6]);
% plotData = struct('points',[],'B',[],'Bu',[],'Bt',[],'Buu',[],...
%     'But',[],'Buut',[]);

% quadData.weightMatrix = spdiags( quadWeightsBT(:), 0, numel(quadWeightsBT), numel(quadWeightsBT));

%quadData.weights = quadWeightsBT(:);
% quadData.weights = {quadWeightsS,quadWeightsT};
% quadData.points = {quadPointsS, quadPointsT};
% %quadData.degree = quadDegree;
% quadData.degree = [nS,nT];
% quadData.knots = knots;
% quadData.noControlPoints = [N,Nt];
%plotData.points = {linspace(0,1,100), linspace(0,1,50) };


if nargout > 1 %Compute outer products for tensor spline
    quadDataTensor = struct('B',[],'Bu',[],'Buu',[],'Bt',[],'But',[],...
        'Buut',[],'quadWeights',[],...
        'Bphi',[],'Bphiu',[]) ;
    quadDataTensor.quadWeights = reshape(quadWeightsS'*quadWeightsT,[],1);
    for ii = N:-1:1
        for jj = Nt:-1:1;
            quadDataTensor.B(:,ii+(jj-1)*N) = reshape(B_S_quad_per(1:4:end,ii)*B_T_quad(1:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1) ;
            quadDataTensor.Bu(:,ii+(jj-1)*N) = reshape(B_S_quad_per(2:4:end,ii)*B_T_quad(1:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Bt(:,ii+(jj-1)*N) = reshape(B_S_quad_per(1:4:end,ii)*B_T_quad(2:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Buu(:,ii+(jj-1)*N) = reshape(B_S_quad_per(3:4:end,ii)*B_T_quad(1:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.But(:,ii+(jj-1)*N) = reshape(B_S_quad_per(2:4:end,ii)*B_T_quad(2:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Buut(:,ii+(jj-1)*N) = reshape(B_S_quad_per(3:4:end,ii)*B_T_quad(2:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            %Third derivative
            quadDataTensor.Buuu(:,ii+(jj-1)*N) = reshape(B_S_quad_per(4:4:end,ii)*B_T_quad(1:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
        end
    end
    
    %phiPath
    for ii = Nphi:-1:1
        for jj = Nt:-1:1;
            
            quadDataTensor.B_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_S_quad_per(1:4:end,ii)*B_T_quad(1:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1) ;
            quadDataTensor.Bu_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_S_quad_per(2:4:end,ii)*B_T_quad(1:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Bt_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_S_quad_per(1:4:end,ii)*B_T_quad(2:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Buu_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_S_quad_per(3:4:end,ii)*B_T_quad(1:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.But_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_S_quad_per(2:4:end,ii)*B_T_quad(2:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Buut_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_S_quad_per(3:4:end,ii)*B_T_quad(2:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Buuu_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_S_quad_per(4:4:end,ii)*B_T_quad(1:2:end,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
        end
    end
    
end

end
