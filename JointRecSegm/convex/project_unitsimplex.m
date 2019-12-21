%
% Author: Jan Lellmann, j.lellmann@damtp.cam.ac.uk
%

function result = project_unitsimplex(c)
%PROJECT_UNITSIMPLEX Projects ROWS of c onto the unit simplex
%   i.e. (euclidean) projection onto the set sum(flatten(x)) = 1, x >= 0.
%   Algorithm converges in at most numel(c) steps, O(n^2) total (could be
%   lower)
%   Source: "C.Michelot: A Finite Algorithm for Finding the Projection of a
%            Point onto the Canonical Simplex of R^n", JOTA Vol.50, No.1,
%            July 1986

    c = c'; % algorithm is formulated for columns

    zeroset = zeros(size(c));
    n = size(c,1);
    ni = n * ones(1,size(c,2));

    x = c;
    
    finished = false;
    while (~finished)
        s = sum(x,1);
        x = (x - repmat((s-1)./ni,[size(c,1) 1])) .* (1-zeroset);
        
        conflicts = find(x < 0);
        if (isempty(conflicts))
            finished = true;
        else
            zeroset(conflicts) = 1;
            ni = n - sum(zeroset,1);
        end
        x(conflicts) = 0;
    end

    result = x'; % algorithm is formulated for columns
end
