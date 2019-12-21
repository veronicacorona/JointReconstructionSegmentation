classdef projection < handle
    
    properties(GetAccess = private, SetAccess = private)        
        dim
        imsize
        sigma
    end
    
    methods(Access = public)
        
        function obj = projection(imsize)
            obj.dim = 2*prod(imsize);
            obj.imsize = [imsize 2];
        end
        
        function sigma = getproxparam(obj)
            sigma = obj.sigma;
        end
        
        function u = prox(obj, f)
            aux = sqrt(sum(abs(reshape(f, [size(obj)/2 2])).^2, 2));
            aux = repmat(aux, [1 2]);
            u = f./max(1, aux(:));
        end
        
        function setproxparam(obj, sigma)
            obj.sigma = sigma;
        end
        
        function dim = size(obj)
            dim = obj.dim;
        end
        
        function imsize = sizend(obj)
            imsize = obj.imsize;
        end
        
    end
    
end