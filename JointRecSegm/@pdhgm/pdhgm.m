classdef pdhgm < handle    
    
    properties(GetAccess = private, SetAccess = private)        
        Fstar
        fig
        G
        K
        k
        maxiter
        ploton = false; 
        res
        resold
        sens
        tol        
        var
        varprev
    end
    
    methods(Access = public)
        
        function obj = pdhgm(K, Fstar, G)
            obj.K = K;
            obj.Fstar = Fstar;
            obj.G = G;
            obj.maxiter = 500;
            obj.tol = 10^-4;
            obj.k = 1;
            obj.sens = 0.001;
        end
        
        %% Additional methods
        
        function disableplot(obj)
            obj.ploton = false;
        end
        
        function enableplot(obj)
            obj.ploton = true;
        end
        
        function initialise(obj)
            obj.k = 1;
            obj.res = Inf;
            obj.resold = 1;
            obj.sens = 0.001;
            obj.var.x = zeros([size(obj.K, 2) 1]);
            obj.var.y = zeros([size(obj.K, 1) 1]);
        end
        
        function restartcounter(obj)
            obj.k = 1;
        end
        
        %% Setter & getter methods
        
        function maxiter = getmaxiter(obj)
            maxiter = obj.maxiter;
        end
        
        function tol = gettolerance(obj)
            tol = obj.tol;
        end
        
          function kit = getiter(obj)
            kit = obj.k;
        end
        
        function var = getvariables(obj)
            var = obj.var;
        end       
        
        function setmaxiter(obj, maxiter)
            obj.maxiter = maxiter;
        end       
        
        function settolerance(obj, tol)
            obj.tol = tol;
        end
        
        function setvariables(obj, var)
            obj.var = var;
        end
        
        %% Main routine
        
        function solve(obj)            
            
            if isempty(obj.var)
                obj.initialise;
            end
            
            if obj.ploton
                obj.fig = figure();              
            end
            
            while ((obj.tol < obj.sens) || (obj.k == 1)) && (obj.k <= ...
                    obj.maxiter)                            
                
                obj.varprev = obj.var;
                
                obj.var.x = obj.G.prox(obj.var.x - obj.G.getproxparam .*...
                    (obj.K'*obj.var.y));
                
                obj.var.y = obj.Fstar.prox(obj.var.y + ...
                     obj.Fstar.getproxparam .* (obj.K*( 2*obj.var.x - ...
                     obj.varprev.x )));
                
                obj.updatesensitivity
                
                if obj.ploton
                    obj.plot
                else
%                     disp(['Iteration nr. ' num2str(obj.k) '\' ...
%                         num2str(obj.maxiter) ', Sensitivity: ' num2str(...
%                         obj.sens)])
                end                
                obj.k = obj.k + 1;                
            end
            
            if obj.ploton
                close(obj.fig)
            else
               % disp('Computation complete!')
            end
            
        end
                
    end
    
    methods(Access = private)        
        
        %% Method for visualisation
                
        function plot(obj)
            obj.fig;
            imagesc(reshape(abs(obj.var.x), sizend(obj.G)))
            axis image
            colormap(gray(512))
            colorbar
            title(['Iteration nr. ' num2str(obj.k) ', Sensitivity: ' ...
                num2str(obj.sens)])
            drawnow
        end
        
        function updatesensitivity(obj)
            aux1 = obj.var.x - obj.varprev.x;
            aux2 = obj.var.y - obj.varprev.y;
            obj.sens = 1/2 * (norm(aux1 - obj.G.getproxparam .* ( ...
                obj.K'*aux2), 2)/norm(obj.var.x, 2) + norm(aux2 - ...
                obj.Fstar.getproxparam .* (obj.K*aux1), 2)/norm( ...
                obj.var.y, 2));
%             obj.sens = 1/2*(norm(obj.var.x - obj.varprev.x, 2)/norm( ...
%                 obj.var.x, 2) + norm(obj.var.y - obj.varprev.y, 2)/norm(...
%                 obj.var.y, 2));
        end
        
    end
    
end