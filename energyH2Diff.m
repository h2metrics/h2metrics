%% energyH2Diff
%
% Computes the energy of dPath after reparametrizing by the linear path
% between phi0 and phi1.
%
% Input
%   dPath
%       Path of curves; matrix of dimension [N*Nt, dSpace]
%   phi0, phi1
%       Reparametrizations of the initial and final curve
%   splineData
%       General information about the splines used.
%   quadData, quadDataTensor
%       Precomputed spline collocation matrices at quadrature points.
%
% Optional input
%   'useDeBoor' = {true (default), false}
%       If true, uses self-written deBoor evaluation algorithm. If false,
%       uses spcol, which is slower, but more reliable.
%
% Output
%   E
%       Energy of the path.
%
function E = energyH2Diff( dPath, phi0, phi1, ...
    splineData, quadData, quadDataTensor, varargin)

useDeBoor = true;

%% Optional parameters
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'usedeboor'
                ii = ii + 1;
                if isnumeric(varargin{ii}) || islogical(varargin{ii})
                    useDeBoor = logical(varargin{ii});
                else
                    error('Invalid value for option ''useDeBoor''.');
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
Nt = splineData.Nt;
Nphi = splineData.Nphi;
nS = splineData.nS;
nT = splineData.nT;
nPhi = splineData.nPhi;

noQuadPointsS = quadData.noQuadPointsS;
noQuadPointsT = quadData.noQuadPointsT;

%% Construct phiPath and evaluate at quadrature sites
phiPath = linearPath(phi0, phi1, splineData);

phiQuad = quadDataTensor.B_phi * phiPath;
phiQuad = reshape(phiQuad, [noQuadPointsS, noQuadPointsT]);
phiQuad = phiQuad + quadData.quadPointsS * ones([1, noQuadPointsT]);

% Check that the columns of phi are indeed diffeomorphisms
phiDiff = diff(phiQuad, 1, 1);
phiDiff = max(phiDiff, splineData.phiEps);
phiQuad = cumsum( [phiDiff(1,:); phiDiff], 1);


%% Evaluate curve at time points
if useDeBoor
    % Use fast method with self-coded deBoor function
    phiQuad = mod(phiQuad, 2*pi);
    cphi_part = zeros([N*noQuadPointsT, dSpace]);
    cphi_t_part = zeros([N*noQuadPointsT, dSpace]);

    B_T = quadData.B_T;
    Bt_T = quadData.Bt_T;

    for kk = N:-1:1
        cphi_part(kk:N:end,:) = B_T * dPath(kk:N:end,:);
        cphi_t_part(kk:N:end,:) = Bt_T * dPath(kk:N:end,:);
    end

    for kk = noQuadPointsT:-1:1
        y = deBoor(splineData.knotsS, nS, ...
            cphi_part((kk-1)*N+1:kk*N,:),phiQuad(:,kk),4,'periodic',1);
        Cphi_u((kk-1)*noQuadPointsS+1:kk*noQuadPointsS,:) = y(2:4:end,:);
        Cphi_uu((kk-1)*noQuadPointsS+1:kk*noQuadPointsS,:) = y(3:4:end,:);
        Cphi_uuu((kk-1)*noQuadPointsS+1:kk*noQuadPointsS,:) = y(4:4:end,:);

        y = deBoor(splineData.knotsS, nS, ...
            cphi_t_part((kk-1)*N+1:kk*N,:),phiQuad(:,kk),3,'periodic',1);
        Cphi_t((kk-1)*noQuadPointsS+1:kk*noQuadPointsS,:) = y(1:3:end,:);
        Cphi_ut((kk-1)*noQuadPointsS+1:kk*noQuadPointsS,:) = y(2:3:end,:);
        Cphi_uut((kk-1)*noQuadPointsS+1:kk*noQuadPointsS,:) = y(3:3:end,:);

    end
else
    % Avoid using deBoor, use spcol instead
    for kk = noQuadPointsT:-1:1
        % Take care of non-periodicity
        phiQuad(:,kk) = mod(phiQuad(:,kk), 2*pi);
        a = find([diff(phiQuad(:,kk)); -1] < 0, 1);
        if a < noQuadPointsS
            if ~isempty(find(diff(phiQuad(a+1:noQuadPointsS,kk), ...
                             1) < 0, 1))
                disp('Should be empty');
            end
            Bphi_S(4*a+1:4*noQuadPointsS,:) = ...
                spcol( splineData.knotsS, nS+1, ...
                    brk2knt( phiQuad(a+1:noQuadPointsS,kk), 4), 'sparse');
        end
        Bphi_S(1:4*a,:) = spcol( splineData.knotsS, nS+1, ...
                                 brk2knt( phiQuad(1:a,kk), 4), 'sparse');

        Bphi_S_per = [Bphi_S(:,1:nS)+Bphi_S(:,end-nS+1:end), ...
            Bphi_S(:,nS+1:end-nS)];
        B_T = quadData.B_T;
        Bt_T = quadData.Bt_T;

        for ii = N:-1:1
            for jj = Nt:-1:1
                Bphi_u((kk-1)*noQuadPointsS+1:kk*noQuadPointsS, ...
                       ii+(jj-1)*N) = Bphi_S_per(2:4:end,ii) * B_T(kk,jj);
                Bphi_uu((kk-1)*noQuadPointsS+1:kk*noQuadPointsS, ...
                       ii+(jj-1)*N) = Bphi_S_per(3:4:end,ii) * B_T(kk,jj);
                Bphi_uuu((kk-1)*noQuadPointsS+1:kk*noQuadPointsS, ...
                       ii+(jj-1)*N) = Bphi_S_per(4:4:end,ii) * B_T(kk,jj);
                Bphi_t((kk-1)*noQuadPointsS+1:kk*noQuadPointsS, ...
                       ii+(jj-1)*N) = Bphi_S_per(1:4:end,ii) * Bt_T(kk,jj);   
                Bphi_ut((kk-1)*noQuadPointsS+1:kk*noQuadPointsS, ...
                       ii+(jj-1)*N) = Bphi_S_per(2:4:end,ii) * Bt_T(kk,jj);
                Bphi_uut((kk-1)*noQuadPointsS+1:kk*noQuadPointsS, ...
                       ii+(jj-1)*N) = Bphi_S_per(3:4:end,ii) * Bt_T(kk,jj);
            end
        end
    end

    Cphi_u = Bphi_u * dPath;
    Cphi_t = Bphi_t * dPath;
    Cphi_ut = Bphi_ut * dPath;
    Cphi_uu = Bphi_uu * dPath;
    Cphi_uut = Bphi_uut * dPath;
    Cphi_uuu = Bphi_uuu * dPath;
end

phi_t = quadDataTensor.Bt_phi * phiPath;
phi_u = 1 + quadDataTensor.Bu_phi * phiPath;
phi_ut = quadDataTensor.But_phi * phiPath;
phi_uu = quadDataTensor.Buut_phi * phiPath;
phi_uut = quadDataTensor.Buut_phi * phiPath;

for ii = dSpace:-1:1
    Cu(:,ii) = Cphi_u(:,ii) .* phi_u;
    Cuu(:,ii) = Cphi_uu(:,ii) .* phi_u .* phi_u + Cphi_u(:,ii) .* phi_uu;
    Ct(:,ii) = Cphi_t(:,ii) + Cphi_u(:,ii) .* phi_t;
    Cut(:,ii) = Cphi_ut(:,ii).*phi_u + Cphi_uu(:,ii).*phi_u.*phi_t + ...
        Cphi_u(:,ii).*phi_ut;
    Cuut(:,ii) = Cphi_uut(:,ii) .* phi_u .* phi_u ...
        + Cphi_uuu(:,ii) .* phi_u .* phi_u .* phi_t ...
        + 2 * Cphi_uu(:,ii) .* phi_u .* phi_ut ...
        + Cphi_ut(:,ii) .* phi_uu + Cphi_uu(:,ii) .* phi_t .* phi_uu ...
        + Cphi_u(:,ii) .* phi_uut;
end

%% Remaining stuff (was identical to energyH2)
Cspeed = sum( Cu.^2 , 2).^(1/2);

% L2 Energy terms
Ct_L2 = sum( Ct.^2,2) .* Cspeed;

% H1 Energy terms
Ct_H1 = sum(Cut.^2, 2) ./ Cspeed;

% H2 Energy terms
CuCuu = zeros([noQuadPointsS*noQuadPointsT, 1]);
CutCut = zeros([noQuadPointsS*noQuadPointsT, 1]);
CutCuut = zeros([noQuadPointsS*noQuadPointsT, 1]);
CuutCuut = zeros([noQuadPointsS*noQuadPointsT, 1]);
for ii = 1:dSpace
    CuCuu = CuCuu + Cu(:,ii) .* Cuu(:,ii);
    CutCut = CutCut + Cut(:,ii) .* Cut(:,ii);
    CutCuut = CutCuut + Cut(:,ii) .* Cuut(:,ii);
    CuutCuut = CuutCuut + Cuut(:,ii) .* Cuut(:,ii);
end

Ct_H2 = CutCut .* CuCuu.^2 ./ Cspeed.^7 ...
    - 2 * CutCuut .* CuCuu ./ Cspeed .^ 5 ...
    + CuutCuut ./ Cspeed .^ 3;

a = splineData.a;
energyIntegrand = a(1) * Ct_L2 + a(2) * Ct_H1 + a(3) * Ct_H2;

%Compute final energy
E = sum(quadDataTensor.quadWeights.*energyIntegrand);

end

