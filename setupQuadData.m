function [quadData, quadDataTensor] = setupQuadData( splineData )
%TODO: Summary
%For optimization
% Input: spData
% Output:
% innerKnotsS, innerKnotsT, quadData struct, (knotS, knotsT)

%clear quadData;
timestamp=now;

quadData = struct('quadPointsS',[],'quadPointsT',[],...
    'noQuadPointsS',[],'noQuadPointsT',[],...,
    'quadWeightsS',[],'quadWeightsT',[],...
    'B_S',[],'Bu_S',[],'Buu_S',[],'Buuu_S',[],...
    'B_T',[],'Bt_T',[],...
    'B_phi',[],'Bu_phi',[],'Buu_phi',[],'Buuu_phi',[],...
    'B_interpolS',[], 'B_interpolPhi',[],'timestamp',timestamp);

quadDegree = splineData.quadDegree;

doS = ~isempty(splineData.nS) && ~isempty(splineData.N) && ...
    ~isempty(splineData.quadDegree);
doT = ~isempty(splineData.nT) && ~isempty(splineData.Nt) && ...
    ~isempty(splineData.quadDegree);
doPhi = ~isempty(splineData.nPhi) && ~isempty(splineData.Nphi);
doInterpolS = ~isempty(splineData.nS) && ~isempty(splineData.N) && ...
    ~isempty(splineData.interpolS);
doInterpolPhi = doPhi && ~isempty(splineData.interpolS);

if doS
    noSubQuadPointsS = ceil((splineData.quadDegree(1)+1)/2);
    N = splineData.N;
    nS = splineData.nS;
    knotsS = splineData.knotsS;
    innerKnotsS = splineData.innerKnotsS;
    
    [quadPointsS, quadWeightsS] = gaussianQuadratureData( ...
        unique(innerKnotsS), 'degree', quadDegree(1) );
    noQuadPointsS = length(quadPointsS);
    
    quadData.quadPointsS = quadPointsS';
    quadData.noQuadPointsS = noQuadPointsS;
    quadData.quadWeightsS = quadWeightsS';
    
    noSder = 3;
    
    B_S_quad = spcol( knotsS, nS+1, ...
                      brk2knt( quadPointsS, noSder+1 ),'sparse');
    B_S_quad_per = [ B_S_quad(:,1:nS) + B_S_quad(:,end-nS+1:end), ...
                     B_S_quad(:,nS+1:end-nS) ];
    
    quadData.B_S = B_S_quad_per(1:4:end,:);
    quadData.Bu_S = B_S_quad_per(2:4:end,:);
    quadData.Buu_S = B_S_quad_per(3:4:end,:);
    quadData.Buuu_S = B_S_quad_per(4:4:end,:);
end

if doT
    noSubQuadPointsT = ceil((splineData.quadDegree(2)+1)/2);
    nT = splineData.nT;
    Nt = splineData.Nt;
    knotsT = splineData.knotsT;
    innerKnotsT = splineData.innerKnotsT;

    [quadPointsT, quadWeightsT] = gaussianQuadratureData( ...
        unique(innerKnotsT), 'degree', quadDegree(2) );

    noQuadPointsT = length(quadPointsT);
    quadData.quadPointsT = quadPointsT';
    quadData.noQuadPointsT = noQuadPointsT;
    quadData.quadWeightsT = quadWeightsT';

    noTder = 1;
    B_T_quad = spcol( knotsT, nT+1, ...
                      brk2knt( quadPointsT, noTder+1 ), 'sparse');              
    quadData.B_T = B_T_quad(1:2:end,:);
    quadData.Bt_T = B_T_quad(2:2:end,:);
end

if doPhi
    nPhi = splineData.nPhi;
    Nphi = splineData.Nphi;
    knotsPhi = splineData.knotsPhi;

    noPhider = 3;
    B_Phi_quad = spcol( knotsPhi, nPhi+1, ...
                        brk2knt( quadPointsS, noPhider+1 ), 'sparse');
    B_Phi_quad_per = [ B_Phi_quad(:,1:nPhi) ...
                        + B_Phi_quad(:,end-nPhi+1:end), ...
                       B_Phi_quad(:,nPhi+1:end-nPhi) ];

    quadData.B_phi = B_Phi_quad_per(1:4:end,:);
    quadData.Bu_phi = B_Phi_quad_per(2:4:end,:);
    quadData.Buu_phi = B_Phi_quad_per(3:4:end,:);
    quadData.Buuu_phi = B_Phi_quad_per(4:4:end,:);
end

if doInterpolS
    N = splineData.N;
    nS = splineData.nS;
    knotsS = splineData.knotsS;
    
    interpolS = splineData.interpolS;
    B_interpolS = spcol( knotsS, nS+1, brk2knt(interpolS, 1), 'sparse');
    B_interpolS = [ B_interpolS(:,1:nS) + B_interpolS(:,end-nS+1:end), ...
                    B_interpolS(:,nS+1:end-nS) ];
    quadData.B_interpolS = B_interpolS;
end

if doInterpolPhi
    B_interpolPhi = spcol( knotsPhi, nPhi+1, ...
                           brk2knt(interpolS, 1), 'sparse');
    B_interpolPhi = [ B_interpolPhi(:,1:nPhi) ...
                        + B_interpolPhi(:,end-nPhi+1:end), ...
                      B_interpolPhi(:,nPhi+1:end-nPhi) ];
    quadData.B_interpolPhi = B_interpolPhi;
end

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

%% Compute outer products for tensor spline
if nargout < 2
    return
end

nnzmaxST = noSubQuadPointsT*(nT+1)*noSubQuadPointsS*(nS+1)*N*Nt;
matrixAllocST = spalloc(noQuadPointsS*noQuadPointsT,N*Nt,nnzmaxST);
%TODO: Are B_phi evaluated at the right points?
%TODO: calculate nnz( quadData.B_phi(:,1)) analytically
nnzmaxPhiT = noSubQuadPointsT*(nT+1)*nnz( quadData.B_phi(:,1) )*Nphi*Nt;
matrixAllocPhiT = spalloc(noQuadPointsS*noQuadPointsT,Nphi*Nt,nnzmaxPhiT);


quadDataTensor = struct('B',matrixAllocST,'Bu',matrixAllocST,'Buu',matrixAllocST,...
    'Bt',matrixAllocST,'But',matrixAllocST,'Buut',matrixAllocST,'Buuu',matrixAllocST,...
    'quadWeights',[],...
    'BuTr',matrixAllocST,'BuuTr',matrixAllocST,'BtTr',matrixAllocST,...
    'ButTr',matrixAllocST,'BuutTr',matrixAllocST,...
    'B_phi',matrixAllocPhiT,'Bu_phi',matrixAllocPhiT,'Bt_phi',matrixAllocPhiT,...
    'Buu_phi',matrixAllocPhiT,'But_phi',matrixAllocPhiT,'Buut_phi',matrixAllocPhiT,...
    'Buuu_phi',matrixAllocPhiT,...
    'timestamp',timestamp,'nnz',[]) ;

% quadDataTensor = struct('B',[],'Bu',[],'Buu',[],'Bt',[],'But',[],...
%     'Buut',[],'quadWeights',[],...
%     'BuTr',[],'BuuTr',[],'BtTr',[],'ButTr',[],'BuutTr',[],...
%     'B_phi',[],'Bu_phi',[],'Bt_phi',[],'Buu_phi',[],'But_phi',[],...
%     'Buut_phi',[],'Buuu_phi',[],'timestamp',timestamp,'nnz',[]) ;

if doS && doT
    quadDataTensor.quadWeights = reshape(quadWeightsS'*quadWeightsT,[],1);
    
%     B_test = spalloc(noQuadPointsS*noQuadPointsT,N*Nt,nnzmaxST);
%     B_test1 = spalloc(noQuadPointsS*noQuadPointsT,N*Nt,nnzmaxST);
    for jj = 1:Nt
        for ii = 1:N
            quadDataTensor.B(:,ii+(jj-1)*N) = reshape(quadData.B_S(:,ii)*quadData.B_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1) ;
            quadDataTensor.Bu(:,ii+(jj-1)*N) = reshape(quadData.Bu_S(:,ii)*quadData.B_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Bt(:,ii+(jj-1)*N) = reshape(quadData.B_S(:,ii)*quadData.Bt_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Buu(:,ii+(jj-1)*N) = reshape(quadData.Buu_S(:,ii)*quadData.B_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.But(:,ii+(jj-1)*N) = reshape(quadData.Bu_S(:,ii)*quadData.Bt_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Buut(:,ii+(jj-1)*N) = reshape(quadData.Buu_S(:,ii)*quadData.Bt_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            %Third derivative
            quadDataTensor.Buuu(:,ii+(jj-1)*N) = reshape(quadData.Buuu_S(:,ii)*quadData.B_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
        end
    end
       
%     for ii = N:-1:1
%         for jj = Nt:-1:1;
% %             B_test1(:,ii+(jj-1)*N) = reshape(B_S_quad_per(1:4:end,ii)*B_T_quad(1:2:end,jj)',...
% %                 noQuadPointsS*noQuadPointsT ,1);
%             
%             quadDataTensor.B(:,ii+(jj-1)*N) = reshape(B_S_quad_per(1:4:end,ii)*B_T_quad(1:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1) ;
%             quadDataTensor.Bu(:,ii+(jj-1)*N) = reshape(B_S_quad_per(2:4:end,ii)*B_T_quad(1:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1)' ;
%             quadDataTensor.Bt(:,ii+(jj-1)*N) = reshape(B_S_quad_per(1:4:end,ii)*B_T_quad(2:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1)' ;
%             quadDataTensor.Buu(:,ii+(jj-1)*N) = reshape(B_S_quad_per(3:4:end,ii)*B_T_quad(1:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1)' ;
%             quadDataTensor.But(:,ii+(jj-1)*N) = reshape(B_S_quad_per(2:4:end,ii)*B_T_quad(2:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1)' ;
%             quadDataTensor.Buut(:,ii+(jj-1)*N) = reshape(B_S_quad_per(3:4:end,ii)*B_T_quad(2:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1)' ;
%             %Third derivative
%             quadDataTensor.Buuu(:,ii+(jj-1)*N) = reshape(B_S_quad_per(4:4:end,ii)*B_T_quad(1:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1)' ;
%         end
%     end

    quadDataTensor.BuTr = quadDataTensor.Bu';
    quadDataTensor.BtTr = quadDataTensor.Bt';
    quadDataTensor.BuuTr = quadDataTensor.Buu';
    quadDataTensor.ButTr = quadDataTensor.But';
    quadDataTensor.BuutTr = quadDataTensor.Buut';
    quadDataTensor.nnz = nnz( quadDataTensor.Bu);
end

%phiPath
if doPhi && doT
    for jj = 1:Nt;
        for ii = 1:Nphi
             quadDataTensor.B_phi(:,ii+(jj-1)*Nphi) = reshape(quadData.B_phi(:,ii)*quadData.B_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1) ;
            quadDataTensor.Bu_phi(:,ii+(jj-1)*Nphi) = reshape(quadData.Bu_phi(:,ii)*quadData.B_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Bt_phi(:,ii+(jj-1)*Nphi) = reshape(quadData.B_phi(:,ii)*quadData.Bt_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Buu_phi(:,ii+(jj-1)*Nphi) = reshape(quadData.Buu_phi(:,ii)*quadData.B_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.But_phi(:,ii+(jj-1)*Nphi) = reshape(quadData.Bu_phi(:,ii)*quadData.Bt_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Buut_phi(:,ii+(jj-1)*Nphi) = reshape(quadData.Buu_phi(:,ii)*quadData.Bt_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
            quadDataTensor.Buuu_phi(:,ii+(jj-1)*Nphi) = reshape(quadData.Buuu_phi(:,ii)*quadData.B_T(:,jj)',...
                noQuadPointsS*noQuadPointsT ,1)' ;
        end
    end
    
%     for ii = Nphi:-1:1
%         for jj = Nt:-1:1;
%             quadDataTensor.B_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_quad_per(1:4:end,ii)*B_T_quad(1:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1) ;
%             quadDataTensor.Bu_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_quad_per(2:4:end,ii)*B_T_quad(1:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1)' ;
%             quadDataTensor.Bt_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_quad_per(1:4:end,ii)*B_T_quad(2:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1)' ;
%             quadDataTensor.Buu_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_quad_per(3:4:end,ii)*B_T_quad(1:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1)' ;
%             quadDataTensor.But_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_quad_per(2:4:end,ii)*B_T_quad(2:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1)' ;
%             quadDataTensor.Buut_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_quad_per(3:4:end,ii)*B_T_quad(2:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1)' ;
%             quadDataTensor.Buuu_phi(:,ii+(jj-1)*Nphi) = reshape(B_Phi_quad_per(4:4:end,ii)*B_T_quad(1:2:end,jj)',...
%                 noQuadPointsS*noQuadPointsT ,1)' ;
%         end
%     end
end
    
end

