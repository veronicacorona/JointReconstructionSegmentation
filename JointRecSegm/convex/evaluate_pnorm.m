%
% Author: Jan Lellmann, j.lellmann@damtp.cam.ac.uk
%

% Takes p-norms row-wise and sums them up.
% This is NOT the weighted lambda/2 TV!!!
% If c is provided, it gives the (possibly negative) weights for the rows.
function [result,resultlocal] = evaluate_pnorm(vecs, p, c)
    %TODO implement infinity-norm
    if (nargin < 3)
        % non-weighted (expanded for performance reasons)
        if (p == 2)
            r = sqrt(sum(vecs.^2,2));
        elseif (p == 1)
            r = sum(abs(vecs),2);
        elseif (p ~= inf)
            r = sum(abs(vecs).^p,2).^(1.0/p);
        else
            error('p = inf not implemented yet');
        end
    else
        % weighted
        if (p == 2)
            r = sqrt(sum(vecs.^2,2)) .* c(:);
        elseif (p == 1)
            r = sum(abs(vecs),2) .* c(:);
        elseif (p ~= inf)
            r = sum(abs(vecs).^p,2).^(1.0/p) .* c(:);
        else
            error('p = inf not implemented yet');
        end
    end
    
    result = sum(r);
    if (nargout > 1)
        resultlocal = r;
    end
end
