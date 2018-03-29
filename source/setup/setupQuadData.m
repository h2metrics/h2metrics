%% setupQuadData
%
% Computes the collocation matrixes for quadrature knots. Computes only
% those matrices, for which information is set in splineData.
%
% Input
%   splineData
%       Information about splines used and quadrature degrees.
%
% Output
%   splineData
%       Adds quadData and quadDataTensor as fields in splineData.
%
function splineData = setupQuadData( splineData )

quadData = struct('quadPointsS', [], 'quadPointsT', [], ...
    'noQuadPointsS', [], 'noQuadPointsT', [], ...
    'quadWeightsS', [], 'quadWeightsT', [], ...
    'B_S', [], 'Bu_S', [], 'Buu_S', [] ,'Buuu_S', [], ...
    'B_T', [], 'Bt_T', [],...
    'B_interpolS', [], 'B_varS', []);

curveClosed = splineData.curveClosed;
quadDegree = splineData.quadDegree;

doS = ~isempty(splineData.nS) && ~isempty(splineData.N) && ...
    ~isempty(splineData.quadDegree);
doT = ~isempty(splineData.nT) && ~isempty(splineData.Nt) && ...
    ~isempty(splineData.quadDegree) && length(splineData.quadDegree) >= 2;
doInterpolS = ~isempty(splineData.nS) && ~isempty(splineData.N) && ...
    ~isempty(splineData.interpolS);
doVar = ~isempty(splineData.varData) && ~isempty(splineData.varData.pts);

if doS
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
    
    noSder = min(4, nS+1);
    
    B_S_quad = spcol( knotsS, nS+1, ...
                      brk2knt( quadPointsS, noSder ), 'sparse');
    if curveClosed
        B_S_quad = [ B_S_quad(:,1:nS) + B_S_quad(:,end-nS+1:end), ...
                     B_S_quad(:,nS+1:end-nS) ];
    end
    
    quadData.B_S = B_S_quad(1:noSder:end, :);
    quadData.Bu_S = B_S_quad(2:noSder:end, :);
    if nS >= 2
        quadData.Buu_S = B_S_quad(3:noSder:end, :);
    end
    if nS >= 3
        quadData.Buuu_S = B_S_quad(4:noSder:end, :);
    end
end

if doT
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

    noTder = 2;
    B_T_quad = spcol( knotsT, nT+1, ...
                      brk2knt( quadPointsT, noTder ), 'sparse');              
    quadData.B_T = B_T_quad(1:noTder:end,:);
    quadData.Bt_T = B_T_quad(2:noTder:end,:);
end



if doInterpolS
    nS = splineData.nS;
    interpolS = splineData.interpolS;
    B_interpolS = spcol( splineData.knotsS, nS+1, ...
                         brk2knt(interpolS, 1), 'sparse');
    if curveClosed
        B_interpolS = [ B_interpolS(:,1:nS) + B_interpolS(:,end-nS+1:end), ...
                        B_interpolS(:,nS+1:end-nS) ];
    end
    quadData.B_interpolS = B_interpolS;
end



%% Compute outer products for tensor spline
quadDataTensor = struct('B',[],'Bu',[],'Buu',[],'Bt',[],'But',[],...
    'Buut',[],'Buuu',[], 'quadWeights',[],...
    'BuTr',[],'BuuTr',[],'BtTr',[],'ButTr',[],'BuutTr',[],...
    'B_phi',[],'Bu_phi',[],'Bt_phi',[],'Buu_phi',[],'But_phi',[],...
    'Buut_phi',[],'Buuu_phi',[],'nnz',[]) ;

if doS && doT
    quadDataTensor.quadWeights = reshape(quadWeightsS'*quadWeightsT,[],1);
    
    quadDataTensor.B = createTensorCollocationMatrix( ...
        quadPointsS, quadPointsT, 1, 1, splineData );
    quadDataTensor.Bu = createTensorCollocationMatrix( ...
        quadPointsS, quadPointsT, 2, 1, splineData );
    quadDataTensor.Bt = createTensorCollocationMatrix( ...
        quadPointsS, quadPointsT, 1, 2, splineData );
    quadDataTensor.But = createTensorCollocationMatrix( ...
        quadPointsS, quadPointsT, 2, 2, splineData );
    
    quadDataTensor.BuTr = quadDataTensor.Bu';
    quadDataTensor.BtTr = quadDataTensor.Bt';
    quadDataTensor.ButTr = quadDataTensor.But';
    quadDataTensor.nnz = nnz( quadDataTensor.Bu);
    
    if nS >= 2
        quadDataTensor.Buu = createTensorCollocationMatrix( ...
            quadPointsS, quadPointsT, 3, 1, splineData );
        quadDataTensor.Buut = createTensorCollocationMatrix( ...
            quadPointsS, quadPointsT, 3, 2, splineData );
        
        quadDataTensor.BuuTr = quadDataTensor.Buu';
        quadDataTensor.BuutTr = quadDataTensor.Buut';
    end
    if nS >= 3
        quadDataTensor.Buuu = createTensorCollocationMatrix( ...
            quadPointsS, quadPointsT, 4, 1, splineData );
    end
end



if doVar
    noPts = splineData.varData.noPts;
    pts = splineData.varData.pts;
    B_varS = spcol( splineData.knotsS, splineData.nS+1, ...
                    brk2knt(pts, 1), 'sparse');
    % Periodicity
    if curveClosed
        nS = splineData.nS;
        B_varS = [B_varS(:,1:nS) + B_varS(:,end-nS+1:end),...
            B_varS(:,nS+1:end-nS) ];
    end
           
    quadData.B_varS = B_varS;
end

%% Save quadData as field in splineData
splineData.quadData = quadData;
splineData.quadDataTensor = quadDataTensor;
    
end