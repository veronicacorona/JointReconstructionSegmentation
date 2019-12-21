classdef dataterm_rec_bregman < handle
    
    properties(GetAccess = public, SetAccess = public)                        
        data
        diagS
        dim
        F        
        tau
        S 
        f
       % C
        pk
        alpha
      %  delta
    end
    
    properties(SetObservable, GetAccess = public, SetAccess = public)                
        proxdata        
    end
    
    methods(Access = public)
        
        function obj = dataterm_rec_bregman(S, F)%, C)
            obj.dim = sizend(F, 2);            
            obj.tau = 1;
            obj.proxdata = zeros([size(S, 1), 1]);
            obj.S = S;
            obj.diagS = diag(S'*S); 
            obj.F = F;
            %obj.C = C;
            obj.pk = zeros(sizend(F, 2));
            obj.alpha = 1;            
            %obj.delta = 1;

            addlistener(obj,'proxdata', 'PostSet', @obj.modifydata);
        end

        function proxdata = getproxdata(obj)
            proxdata = obj.proxdata;
        end
        
        function tau = getproxparam(obj)
            tau = obj.tau;
        end        
        function alpha = getproxparam1(obj)
            alpha = obj.alpha;
        end
        function delta = getproxparam2(obj)
            delta = obj.delta;
        end
        function modifydata(obj, ~, ~)
            obj.data = obj.S'*obj.proxdata;
        end
        
        function u = prox(obj, f)
            
            %u = obj.F'*((obj.F*f(:) + obj.getproxparam .* obj.data)./(1 ...
             %   + obj.getproxparam .* obj.diagS));
         u = obj.F'*((obj.F*(f(:)+obj.getproxparam.*obj.getproxparam1.*obj.pk) + obj.getproxparam .* obj.data)./(1 ...
                + obj.getproxparam .* obj.diagS));
        
        end
        
        function setproxdata(obj, proxdata)
            obj.proxdata = proxdata(:);
        end
         function setC(obj, C)
            obj.C = C(:);
         end
         function C = getC(obj)
             C= obj.C;
            
         end
         function setP(obj, pk)
            obj.pk = pk(:);
         end
         function pk = getP(obj)
             pk= obj.pk;
            
         end
        function setproxparam(obj, tau)
            obj.tau = tau;            
        end
         function setproxparam1(obj, alpha)
            obj.alpha = alpha;            
         end
        function setproxparam2(obj, delta)
            obj.delta = delta;            
        end
        function dim = size(obj)
            dim = prod(obj.dim);
        end
        
        function dim = sizend(obj)
            dim = obj.dim;
        end
        
    end
    
end