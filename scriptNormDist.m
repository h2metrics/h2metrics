%Load results and compute norm distances
%Change these values as needed

load('spline_results.mat');

nS = 4;
nT  = 1;
quadDegree = [8,4];

%Pairs of (N,Nt) to test
test_list = [ (10:10:70)' ,repmat( 12, 7,1) ]; 

task = 'circle2transcircle';
% Add more task names

normMat = ones( 2, length(Nt_list)-1 );

for i = 1:size(test_list,1)-1
        disp(['i = ',num2str(i)]);
        %First result
        resultname1 = [task,'_','n',num2str(nS),'_','nT',num2str(nT),...
            '_','N',num2str(test_list(i,1)),...
            '_','Nt',num2str(test_list(i,2)),...
            'q',num2str(quadDegree(1)),num2str(quadDegree(2))];
        %Second result, j + 1
        resultname2 = [task,'_','n',num2str(nS),'_','nT',num2str(nT),...
            '_','N',num2str(test_list(i+1,1)),...
            '_','Nt',num2str(test_list(i+1,2)),...
            'q',num2str(quadDegree(1)),num2str(quadDegree(2))];
        
        [H1H2norm, L2L2norm] = pathNormDist(spline_results.(resultname1),...
            spline_results.(resultname2) );
        
        normMat( 1, i) = H1H2norm;
        normMat( 2, i) = L2L2norm;
        
end

semilogy( normMat' )
ylabel('||c_i - c_{(i-1)}||')
legend('H1H2','L2L2')
title('Norm difference between solutions')

%% Plot the pointwise energy of a path

nS = 4;
nT  = 1;
N = 70;
Nt = 12;

resultname = [task,'_','n',num2str(nS),'_','nT',num2str(nT),...
    '_','N',num2str(N),...
    '_','Nt',num2str(Nt),...
    'q',num2str(quadDegree(1)),num2str(quadDegree(2))];

tempstruct = spline_results.(resultname);
tempstruct.('Nt_inner') = tempstruct.Nt - 2;

tempstruct.('knotsS') =  [(-tempstruct.nS):(tempstruct.N+tempstruct.nS)]/tempstruct.N*2*pi; %normalize, domain of definition is [0,2*pi]x[0,2*pi]
tempstruct.('knotsT') = [ zeros(1,tempstruct.nT), linspace(0,1,tempstruct.Nt - tempstruct.nT + 1), ones(1,tempstruct.nT)];

tempstruct.('quadDegree') = [4,4];

[innerKnotsS, innerKnotsT, quadData] = setupQuadData(tempstruct);
speed = evaluateSpeed(tempstruct.d, quadData);

plot( linspace(0,1,length(speed)), speed)
