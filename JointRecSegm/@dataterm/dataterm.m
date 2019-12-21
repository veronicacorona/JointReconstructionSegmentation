
classdef dataterm < handle
    
    properties(GetAccess = private, SetAccess = private)                        
        data
        diagS
        dim
        F        
        tau
        S                
    end
    
    properties(SetObservable, GetAccess = private, SetAccess = private)                
        proxdata        
    end
    
    methods(Access = public)
        
        function obj = dataterm(S, F)
            obj.dim = sizend(F, 2);            
            obj.tau = 1;
            obj.proxdata = zeros([size(S, 1), 1]);
            obj.S = S;
            obj.diagS = diag(S'*S); 
            obj.F = F;
            addlistener(obj,'proxdata', 'PostSet', @obj.modifydata);
        end

        function proxdata = getproxdata(obj)
            proxdata = obj.proxdata;
        end
        
        function tau = getproxparam(obj)
            tau = obj.tau;
        end        
        
        function modifydata(obj, ~, ~)
            obj.data = obj.S'*obj.proxdata;
        end
        
        function u = prox(obj, f)
            u = obj.F'*((obj.F*f(:) + obj.getproxparam .* obj.data)./(1 ...
                + obj.getproxparam .* obj.diagS));
        end
        
        function setproxdata(obj, proxdata)
            obj.proxdata = proxdata(:);
        end
        
        function setproxparam(obj, tau)
            obj.tau = tau;            
        end
        
        function dim = size(obj)
            dim = prod(obj.dim);
        end
        
        function dim = sizend(obj)
            dim = obj.dim;
        end
        
    end
    
end