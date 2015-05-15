%% energyH2Diff
%
% Computes the energy of dPath after some reparametrizations...
%
% Input
%   dPath
%       Path of curves; matrix of dimension [N*Nt, dSpace]
%   phi
%       Reparametrization of the final curve
%   splineData
%       General information about the splines used.
%   quadData, quadDataTensor
%       Precomputed spline collocation matrices at quadrature points.
%
% Output
%   E
%       Energy of the path.
%
function E = energyH2Diff( dPath, phi, v, beta, ...
    splineData, quadData, quadDataTensor, varargin)

optDiff = true;
optTra = true;
optRot = true;

%% Optional parameters
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'optdiff'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    optDiff = logical(varargin{ii});
                else
                    error('Invalid value for option ''optDiff''.');
                end
            case 'opttra'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    optTra = logical(varargin{ii});
                else
                    error('Invalid value for option ''optTra''.');
                end
            case 'optrot'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    optRot = logical(varargin{ii});
                else
                    error('Invalid value for option ''optRot''.');
                end
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1; 
    end
end

%% Extract parameters
dSpace = splineData.dSpace;
N = splineData.N;

%% Apply diffeomorphism, translation and rotation
d1 = dPath(end-N+1:end,:);

if optDiff
    d1 = curveComposeDiff(d1, phi, splineData, quadData);
end
if optTra
    d1 = d1 + ones([N, 1]) * v';
end
if optRot
    rotation = [ cos(beta), -sin(beta); ...
                 sin(beta),  cos(beta) ];
    d1 = d1 * rotation;
end

dPath2 = dPath;
dPath2(end-N+1:end,:) = d1;

%% Evaluate path at quadrature sites
Cu = quadDataTensor.Bu * dPath2;
Ct = quadDataTensor.Bt * dPath2;
Cut = quadDataTensor.But * dPath2;
Cuu = quadDataTensor.Buu * dPath2;
Cuut = quadDataTensor.Buut * dPath2;

%% Remaining stuff (should be identical to energyH2)
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

a = splineData.a;
energyIntegrand = a(1) * Ct_L2 + a(2) * Ct_H1 + a(3) * Ct_H2;

%Compute final energy
E = sum(quadDataTensor.quadWeights.*energyIntegrand);

end

