function analyzeBvpResults( d1, d2, dPath, splineData, options, info )

disp(datestr(now, 'HH:MM:SS'));

%% Basic parameters
disp(' ');
disp('Basic parameters');

if isfield(splineData, 'scaleInv')
    scaleInv = splineData.scaleInv;
else
    scaleInv = 0;
end
disp(['scaleInv =', num2str(scaleInv)]);
disp(['a =', num2str(splineData.a)]);
    
%% Curve length
disp(' ');
disp(['Curve 1, Len=', ...
      num2str(curveLength(d1, splineData))]);
disp(['Curve 2, Len=', ...
      num2str(curveLength(d2, splineData))]);
  
%% Geodesic distance
en12 = pathRiemH2Energy(dPath, splineData);
dist12 = sqrt(en12);
varEnd2 = varifoldDistanceSquared(...
    evalPath(dPath, 1, splineData), d2, splineData);

disp(' ');
disp('Geodesic distance and varifolds');
disp(['dist(d1,d2)   =', num2str(dist12)]);
disp(['dist(d1,d2)^2 =', num2str(en12)]);
disp(['var(d(1),d2)^2=', num2str(varEnd2)]);
disp(['lambda * var(d(1),d2)^2=', num2str(options.varLambda * varEnd2)]);
disp(['var(d1,d2)^2  =', ...
      num2str(varifoldDistanceSquared(d1, d2, splineData))]);
disp(['dist(d1,d2)^2 / var(d(1),d2)^2 = ', num2str(en12 / varEnd2)]);
disp(['lambda =', num2str(options.varLambda)]);


%% The parametrized problem
% dPathEnd = evalPath(dPath, 1, splineData);
% optParam = options;
% optParam.optDiff = false;
% optParam.useMultigrid = false;
% optParam.hansoPrtLevel = 0;
% [~, optPathParam, ~, ~] = geodesicBvp(d1, dPathEnd, splineData, ...
%                                           optParam, 'initPath', dPath);
% 
% enParam = pathRiemH2Energy(optPathParam, splineData);
% distParam = sqrt(enParam);
% 
% disp(' ');
% disp('Comparing d1, d(1) with parametrized problem.');
% disp(['dist(d1,dp(1))   =', num2str(distParam)]);
% disp(['dist(d1,dp(1))^2 =', num2str(enParam)]);

%% Split metric into components
disp(' ');
disp('Components of the path energy');
[E, comp] = pathRiemH2Energy(dPath, splineData);
comp = comp / E;
disp(['Rel. weight =', num2str(comp, 2)]);

%% Check that path does not have self-intersections
[cpMin, cpMax, turningSame] = checkPathImmersion(dPath, splineData);

disp(' ');
disp(['Along geodesic: min |c''| = ', num2str(cpMin)]);
disp(['Along geodesic: max |c''| = ', num2str(cpMax)]);
disp(['Turing number remains same: ', num2str(turningSame)]);

%% Info about optimization
if isempty(info) || ~isfield(info, 'infoHanso')
    return
end
    
infoHanso = info.infoHanso;

disp(' ');
disp('Information from Hanso');
disp(['iter =', num2str(infoHanso.iter)]);
disp(['dnorm =', num2str(infoHanso.dnorm)]);
disp(['riemH1H2norm of grad =', ...
    num2str(pathRiemH1H2InnerProd(dPath, ...
        infoHanso.grad, infoHanso.grad, splineData))]);
disp(['duration =', num2str(infoHanso.cputime)]);
disp(['stopping criteria =', num2str(infoHanso.info)]);

end

