%
% Author: Jan Lellmann, j.lellmann@damtp.cam.ac.uk
%

% EUCLIDEAN projection of all rows (!) in input to unit sphere of radius 1en in p-norm
% currently implemented for p=2, p=inf
% WARNING: for p other than 2,inf, projection is NOT trivial!
function result = project_pnorm_uniform(vectors, len, p)
  if (p == 2)
      norms = sqrt(sum(vectors.^2,2)); % row vector with norms
      factors = 1./max(1,norms/len);  
      result = spdiag(factors) * vectors; %SLOW
  elseif (p == inf)
      norms = abs(vectors);
      result = vectors ./ max(1,norms/len); %SLOW
  else
      error('p must be 2 or +inf');
  end
end
