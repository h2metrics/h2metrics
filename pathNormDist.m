function [ H1H2norm, L2L2norm ] = pathNormDist( splinePath1, splinePath2)
%Calculate different measures of distance between two spline paths.

%Compute knot sequence in space and time
knotsS1 =  [(-splinePath1.nS):(splinePath1.N+splinePath1.nS)]/splinePath1.N*2*pi; %normalize, domain of definition is [0,2*pi]x[0,2*pi]
knotsT1 = [ zeros(1,splinePath1.nT), linspace(0,1,splinePath1.Nt - splinePath1.nT + 1), ones(1,splinePath1.nT)];
innerKnotsS1 = knotsS1(splinePath1.nS+1:end-splinePath1.nS);
innerKnotsT1 = knotsT1(splinePath1.nT+1:end-splinePath1.nT);

knotsS2 =  [(-splinePath2.nS):(splinePath2.N+splinePath2.nS)]/splinePath2.N*2*pi; %normalize, domain of definition is [0,2*pi]x[0,2*pi]
knotsT2 = [ zeros(1,splinePath2.nT), linspace(0,1,splinePath2.Nt - splinePath2.nT + 1), ones(1,splinePath2.nT)];
innerKnotsS2 = knotsS2(splinePath2.nS+1:end-splinePath2.nS);
innerKnotsT2 = knotsT2(splinePath2.nT+1:end-splinePath2.nT);

%Find intervals where both piecewise smooth on
knotsS = sort(unique([innerKnotsS1 innerKnotsS2]));
knotsT = sort(unique([innerKnotsT1 innerKnotsT2]));

%Quadrature sites and weights
[quadPointsS, weightsS] = gaussianQuadratureData( knotsS, 'degree',4 );
[quadPointsT, weightsT] = gaussianQuadratureData( knotsT, 'degree',4 );
noQuadPointsS = length(quadPointsS);
noQuadPointsT = length(quadPointsT);

%Evaluate basis functions, periodic in space
B_S1 = spcol( knotsS1, splinePath1.nS+1, brk2knt( quadPointsS, 3 ),'sparse');
B_T1 = spcol( knotsT1, splinePath1.nT+1, brk2knt( quadPointsT, 2 ),'sparse');
B_S1_per = [B_S1(:,1:splinePath1.nS) + B_S1(:,end-splinePath1.nS+1:end),...
    B_S1(:,splinePath1.nS+1:end-splinePath1.nS)];

B_S2 = spcol( knotsS2, splinePath2.nS+1, brk2knt( quadPointsS, 3 ),'sparse');
B_T2 = spcol( knotsT2, splinePath2.nT+1, brk2knt( quadPointsT, 2 ),'sparse');
B_S2_per = [B_S2(:,1:splinePath2.nS) + B_S2(:,end-splinePath2.nS+1:end),...
    B_S2(:,splinePath2.nS+1:end-splinePath2.nS)];

%Evaluate spline paths
for ii = splinePath1.N:-1:1
    for jj = splinePath1.Nt:-1:1;
    B1(:,ii+(jj-1)*splinePath1.N) = reshape(B_S1_per(1:3:end,ii)*B_T1(1:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    B1s(:,ii+(jj-1)*splinePath1.N) = reshape(B_S1_per(2:3:end,ii)*B_T1(1:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    B1ss(:,ii+(jj-1)*splinePath1.N) = reshape(B_S1_per(3:3:end,ii)*B_T1(1:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    B1t(:,ii+(jj-1)*splinePath1.N) = reshape(B_S1_per(1:3:end,ii)*B_T1(2:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    B1st(:,ii+(jj-1)*splinePath1.N) = reshape(B_S1_per(2:3:end,ii)*B_T1(2:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    B1sst(:,ii+(jj-1)*splinePath1.N) = reshape(B_S1_per(3:3:end,ii)*B_T1(2:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    end
end

C1 = B1*splinePath1.d;
C1s = B1s*splinePath1.d;
C1ss = B1ss*splinePath1.d;
C1t = B1t*splinePath1.d;
C1st = B1st*splinePath1.d;
C1sst = B1sst*splinePath1.d;

for ii = splinePath2.N:-1:1
    for jj = splinePath2.Nt:-1:1;
    B2(:,ii+(jj-1)*splinePath2.N) = reshape(B_S2_per(1:3:end,ii)*B_T2(1:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    B2s(:,ii+(jj-1)*splinePath2.N) = reshape(B_S2_per(2:3:end,ii)*B_T2(1:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    B2ss(:,ii+(jj-1)*splinePath2.N) = reshape(B_S2_per(3:3:end,ii)*B_T2(1:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    B2t(:,ii+(jj-1)*splinePath2.N) = reshape(B_S2_per(1:3:end,ii)*B_T2(2:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    B2st(:,ii+(jj-1)*splinePath2.N) = reshape(B_S2_per(2:3:end,ii)*B_T2(2:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    B2sst(:,ii+(jj-1)*splinePath2.N) = reshape(B_S2_per(3:3:end,ii)*B_T2(2:2:end,jj)',...
        noQuadPointsS*noQuadPointsT ,1) ;
    end
end

C2 = B2*splinePath2.d;
C2s = B2s*splinePath2.d;
C2ss = B2ss*splinePath2.d;
C2t = B2t*splinePath2.d;
C2st = B2st*splinePath2.d;
C2sst = B2sst*splinePath2.d;

%Construct weight matrix
weightMult = (weightsT'*weightsS)'; %Multiply u t(i)
weightMatrix = spdiags( weightMult(:), 0, numel(weightMult), numel(weightMult));

%Evaluate norms
%H1H2
H1H2norm = sum(weightMatrix*( sum( (C1-C2).^2,2) + sum( (C1s-C2s).^2,2) + ...
    sum( (C1ss-C2ss).^2,2) + sum( (C1t-C2t).^2,2) + sum( (C1st-C2st).^2,2) + ...
    sum( (C1sst-C2sst).^2,2) ));

L2L2norm = sum( weightMatrix*( sum( (C1-C2).^2,2) + sum( (C1t-C2t).^2,2) ) );

end

