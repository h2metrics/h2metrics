function [E_geo, dPath, opt_output, opt_grad] = geodesicBVP(d0,d1,splineData,quadData,quadDataTensor,varargin);
% Compute the minimal geodesic connecting the splines given by d0 and d1.%
%
% Input: 
%       d0, [NxdSpace], first set of control points
%       d1, [NxdSpace], second set of control points
%       splineData,
%       quadData,
%       quadDataTensor
%       'option' = 'value', list of options on the form 
%       
% Output:
%       E_geo, optimal energy
%       dPath, optimal path
%
options = optimoptions('fminunc');
options = optimoptions(options, 'MaxFunEvals',300000);

d_linear = linearPath(d0,d1,splineData.Nt);
d_init = d_linear( splineData.N+1:end-splineData.N,:);

%Some code for handling optional inputs
ii = 1;
while ii <= length(varargin)-1
    if (isa(varargin{ii},'char') && isa(varargin{ii+1},'char'))
         switch (lower(varargin{ii}))
             case 'display' %Don't return gradient terms from end curves
                 options = optimoptions(options,'Display',lower(varargin{ii+1}));
             case 'plotfval' %value = 'true'/'false'
                 if strcmpi(varargin{ii+1},'true')
                    options = optimoptions(options,'PlotFcns', @optimplotfval); 
                 end
             case 'algorithm'
                 options = optimoptions(options,'Algorithm',lower(varargin{ii+1}));
             case 'gradobj'
                 options = optimoptions(options,'GradObj', lower(varargin{ii+1}));
         end
    elseif (isa(varargin{ii},'char') && isa(varargin{ii+1},'numeric'))
        switch (lower(varargin{ii}))
            case 'tolfun'
                options = optimoptions(options,'TolFun', varargin{ii+1});
            case 'tolx'
                options = optimoptions(options,'TolX', varargin{ii+1});
            case 'init'
                d_init = varargin{ii+1};
            case 'MaxFunEvals'
                options = optimoptions(options, 'MaxFunEvals',varargin{ii+1});
        end
    end
    ii = ii + 1;  
end

[d_optimal,E_geo,exitflag,output,grad] = ...
fminunc(@(d_opt)energyH2([d0;d_opt;d1],splineData,quadData,quadDataTensor),...
    d_init,options);

dPath = [d0;d_optimal;d1];

if nargout > 2 %provide output options
    opt_output = output;
end
if nargout > 3 %provide output options
    opt_grad = grad;
end

%% TODO: Some hessian options

%options = optimoptions(options,'Hessian', 'off');
%options = optimoptions(options,'Algorithm', 'quasi-newton');
%options = optimoptions(options,'HessUpdate', 'dfp');

% Uncomment three next lines to use trust-region with simplified Hessian
% options = optimoptions(options,'Algorithm', 'trust-region');
% hessPattern = HessianPattern(nS,nT,N,Nt);
% options = optimoptions(options,'HessPattern', hessPattern);

%return hessian



end

