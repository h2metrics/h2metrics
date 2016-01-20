function [G, comp] = pathRiemH2Length( dPath, splineData, varargin )
                                   
% Set constants to be used in metric
a = [1 0 1];
if ~isempty(splineData.a)
    a = splineData.a;
end

ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'a'
                ii = ii + 1;
                a = varargin{ii};
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1; 
    end
end

dSpace = splineData.dSpace;
quadData = splineData.quadData;
quadDataTensor = splineData.quadDataTensor;

%% Evaluate path at quadrature sites
Cu = quadDataTensor.Bu * dPath;
Ct = quadDataTensor.Bt * dPath;
Cut = quadDataTensor.But * dPath;
Cuu = quadDataTensor.Buu * dPath;
Cuut = quadDataTensor.Buut * dPath;

%% Calculate terms of the energy
Cspeed = Cu(:,1) .* Cu(:,1);
for ii = 2:dSpace
    Cspeed = Cspeed + Cu(:,ii) .* Cu(:,ii);
end
Cspeed = sqrt(Cspeed);
% Cspeed = sqrt( sum( Cu.^2 , 2) );
CspeedInv = 1 ./ Cspeed;

% L2 and H1 terms
Ct_L2 = Ct(:,1) .* Ct(:,1);
Ct_H1 = Cut(:,1) .* Cut(:,1);
for ii = 2:dSpace
    Ct_L2 = Ct_L2 + Ct(:,ii) .* Ct(:,ii);
    Ct_H1 = Ct_H1 + Cut(:,ii) .* Cut(:,ii);
end
Ct_L2 = Ct_L2 .* Cspeed;
Ct_H1 = Ct_H1 .* CspeedInv;

% Ct_L2 = sum( Ct .* Ct, 2) .* Cspeed;
% Ct_H1 = sum( Cut .* Cut, 2) .* CspeedInv;

% H2 Energy terms
CuCuu = Cu(:,1) .* Cuu(:,1);
CutCut = Cut(:,1) .* Cut(:,1);   
CutCuut = Cut(:,1) .* Cuut(:,1);
CuutCuut = Cuut(:,1) .* Cuut(:,1);

for ii = 2:dSpace
    CuCuu = CuCuu + Cu(:,ii) .* Cuu(:,ii);
    CutCut = CutCut + Cut(:,ii) .* Cut(:,ii);
    CutCuut = CutCuut + Cut(:,ii) .* Cuut(:,ii);
    CuutCuut = CuutCuut + Cuut(:,ii) .* Cuut(:,ii);
end

% Ct_H2 = CutCut .* CuCuu.^2 ./ Cspeed.^7 ...
%     - 2 * CutCuut .* CuCuu ./ Cspeed .^ 5 ...
%     + CuutCuut ./ Cspeed .^ 3;

CspeedInv2 = CspeedInv .* CspeedInv;
CspeedInv3 = CspeedInv2 .* CspeedInv;
CspeedInv5 = CspeedInv3 .* CspeedInv2;
CspeedInv7 = CspeedInv5 .* CspeedInv2;
Ct_H2 = CutCut .* CuCuu.^2 .* CspeedInv7 ...
    - 2 * CutCuut .* CuCuu .* CspeedInv5 ...
    + CuutCuut .* CspeedInv3;

%% We want to compute length, not energy
noQuadPointsS = quadData.noQuadPointsS;
noQuadPointsT = quadData.noQuadPointsT;
for jj = noQuadPointsT:-1:1
    Ct_L2_T(jj) = sum( Ct_L2((jj-1)*noQuadPointsS+1:jj*noQuadPointsS) ...
        .* quadData.quadWeightsS );
    Ct_H1_T(jj) = sum( Ct_H1((jj-1)*noQuadPointsS+1:jj*noQuadPointsS) ...
        .* quadData.quadWeightsS );
    Ct_H2_T(jj) = sum( Ct_H2((jj-1)*noQuadPointsS+1:jj*noQuadPointsS) ...
        .* quadData.quadWeightsS );
    
    Ct_L2_T(jj) = sqrt(Ct_L2_T(jj));
    Ct_H1_T(jj) = sqrt(Ct_H1_T(jj));
    Ct_H2_T(jj) = sqrt(Ct_H2_T(jj));
    dist_T(jj) = sqrt(a(1) * Ct_L2_T(jj) + a(2) * Ct_H1_T(jj) + ...
        a(3) * Ct_H2_T(jj));
end
dist_L2 = sum(Ct_L2_T' .* quadData.quadWeightsT);
dist_H1 = sum(Ct_H1_T' .* quadData.quadWeightsT);
dist_H2 = sum(Ct_H2_T' .* quadData.quadWeightsT);
G = sum(dist_T' .* quadData.quadWeightsT);

if nargout > 1
    comp = [sqrt(a(1))*dist_L2, sqrt(a(2))*dist_H1, sqrt(a(3))*dist_H2];
end