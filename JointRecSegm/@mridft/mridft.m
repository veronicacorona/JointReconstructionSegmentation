classdef mridft < fcthdlop          
    
    methods(Access = public)
        function obj = mridft(domaindim, center)                        
            if (nargin == 0)
                domaindim = [256 256];                
            end
            if (nargin <= 1)
                if (numel(domaindim) == 2) && (domaindim(2) == 1)
                    center = ceil((domaindim(1) + 1)/2);
                else                    
                    center = ceil((domaindim + ones(size(domaindim))/2));     
                end
            end
            if ~isequal(numel(domaindim), numel(center))
                error('Dimensions of signal and centerpoint do not match!')
            end
            if numel(center) == 1
                onevec = 1;
            else
                onevec = ones(size(domaindim));
            end
            fwfcthdl = @(u) sqrt(numel(u))*circshift((ifftn(ifftshift(u...
                ))), center - onevec);
            bwfcthdl = @(f) 1/sqrt(numel(f))*fftshift(fftn(circshift(f, ...
                -center + onevec)));
            obj = obj@fcthdlop(domaindim, domaindim, fwfcthdl, bwfcthdl);            
        end
    end       
    
end