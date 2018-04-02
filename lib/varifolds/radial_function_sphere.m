function [res, dres]=radial_function_sphere(x,objfun)
% r=RHO_SPH(x,opt,sig) implements the Gaussian kernel and its derivatives :
%
%       rho(t)=exp(-t/lam^2);
%
% Computation of rho(x^2/2) (ie opt==0), rho'(x^2/2) (ie opt==1) and
% rho"(x^2/2)  (ie opt==2)
%
% Input :
%   x : a matrix
%   opt : a integer
%   sig : vector with kernel bandwidths
%
% Output :
%   r : matrix of the size of x
%
% See also : .
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, G. Nardi, A. Trouve (2012-2016)

res = zeros(size(x));
dres = zeros(size(x));

switch lower(objfun.kernel_grass)
    case 'linear'
        res = x;
        dres = ones(size(x));
        
%         if (derivative_order==0)
%             res=x;
%         elseif (derivative_order==1)
%             res=ones(size(x));
%         end
        
    case 'gaussian_oriented'
        
        for l = objfun.kernel_size_grass
            e1 = exp(-2/l^2) * exp(2*x/l^2);
            res = res + e1;
            dres = dres + 2 * e1 / l^2;
            
%             if derivative_order==0
%                 res = res + exp(-2/l^2) * exp(2*x/l^2);
%             elseif derivative_order==1
%                 res = res + 2 * exp(-2/l^2) * exp(2*x/l^2) / l^2;
%             end
        end
        
        
    case 'binet'
        
        res = x.^2;
        dres = 2*x;
        
%         if (derivative_order==0)
%             res=x.^2;
%         elseif (derivative_order==1)
%             res=2*x;
%         end
        
    case 'gaussian_unoriented'
        
        e1 = exp(-2/l^2) * exp(2*x.^2/l^2);
        res = res + e1;
        dres = dres + 4 * x .* e1 / l^2;
        
%         for l=objfun.kernel_size_grass
%             if (derivative_order==0)
%                 res = res + exp(-2/l^2) * exp(2*x.^2/l^2);
%             elseif (derivative_order==1)
%                 res = res + (4 * exp(-2/l^2) / l^2) * x .* exp(2*x.^2/l^2);
%             end
%         end
        
    otherwise
        
        res = objfun.distance_kernel(x);
        dres = objfun.distance_kernel_der(x);
        
%         if (derivative_order==0)
%             res=objfun.distance_kernel(x);
%         elseif (derivative_order==1)
%             res=objfun.distance_kernel_der(x);
%         end
        
end
end

