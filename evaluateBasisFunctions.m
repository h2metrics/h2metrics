function [ B ] = evaluateBasisFunctions( evalS,evalT,spData, varargin )
%Evaluate the Basis functions in the tensor product B-spline given by
%spData, at the sites given in the vectors s and t.

Der = [0,0];
i = 1;
while i <= length(varargin)
    if (isa(varargin{i},'char'))
         switch (lower(varargin{i}))
             case 'der' %Compute derivatives
                 Der = varargin{i+1};
             otherwise
                 error('Invalid option: ''%s''.',varargin{i});
         end   
    end
    i = i + 1;  
end

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

B_S = spcol( knotsS, nS+1, brk2knt( evalS, Der(1)+1 ),'sparse');
B_T = spcol( knotsT, nT+1, brk2knt( evalT, Der(2)+1 ),'sparse');

B_S_per = [B_S(:,1:nS) + B_S(:,end-nS+1:end), B_S(:,nS+1:end-nS)];

for ii = N:-1:1
    for jj = Nt_inner+2:-1:1;
    B(:,ii+(jj-1)*N) = reshape(B_S_per(1:Der(1)+1:end,ii)*B_T(1:Der(2)+1:end,jj)',...
        length(evalS)*length(evalT) ,1) ;
    
    end
end

end

