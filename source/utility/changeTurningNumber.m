function change = changeTurningNumber(dPath, splineData, varargin)

% Handle optional inputs
p = inputParser;
addParameter(p, 'noPtsT', splineData.Nt);
parse(p, varargin{:});

noPtsT = p.Results.noPtsT;
ptsT = linspace(0, 1, noPtsT);

quadData = splineData.quadData;

%% Calculate the turning number
% total curvature is total integral of curvature / 2*pi
% D^2 c = D (c'/|c'|) = c'' / |c'|^2 - c' <c'',c'> / |c'|^4
% <D^2 c,n> ds = <D^2 c, Jc'> = <c'', Jc'> / |c'|^2
% D^2 c = kn

turningNumber = zeros(noPtsT, 1);
for jj = 1:noPtsT
    d = evalPath(dPath, ptsT(jj), splineData);
    
    Cu = quadData.Bu_S*d;
    Cuu = quadData.Buu_S*d;

    Cspeed = sum( Cu.^2 , 2).^(1/2);
    CspeedInv = 1./Cspeed;

    J = [ 0, 1; ...
         -1, 0 ];
    JCu = Cu * J;
    CuuJCu = sum(Cuu .* JCu, 2);

    curvatureIntegrand = CuuJCu .* CspeedInv.^2;
    turningNumber(jj) = ...
        quadData.quadWeightsS' * curvatureIntegrand / (2*pi);
    turningNumber(jj) = round(turningNumber(jj));
end

change = ~all(turningNumber == turningNumber(1));

end