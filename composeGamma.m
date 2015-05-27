function ga = composeGamma(ga1, ga2, splineData, quadData)

ga = struct( 'phi', [], 'beta', [], 'v', [], 'alpha', []);

ga.v = ga1.v;
if ~isempty(ga1.beta) && ~isempty(ga2.v)
    rotation = [ cos(ga1.beta), -sin(ga1.beta) ...
                 sin(ga1.beta), cos(ga1.beta) ];
	ga.v = ga.v + rotation * ga2.v;
elseif ~isempty(ga2.beta)
    ga.v = ga.v + ga2.b;

ga.beta = ga1.beta;
if ~isempty(ga2.beta)
    ga.beta = ga.beta + ga2.beta;
end

ga.alpha = ga1.alpha;
if isempty(ga1.phi) && ~isempty(ga2.alpha)
    ga.alpha = ga.alpha + ga2.alpha;
end

ga.phi = ga1.phi;
if ~isempty(ga2.phi)
    if ~isempty(ga2.alpha)
        ga.phi = ga2.phi + ...
            curveComposeDiff( ga1.phi, ga2.phi - ga2.alpha, ...
                              splineData, quadData );
    else
        ga.phi = ga2.phi + ...
            curveComposeDiff( ga1.phi, ga2.phi, splineData, quadData );
    end
end


    