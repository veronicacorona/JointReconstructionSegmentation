%
% Author: Jan Lellmann, j.lellmann@damtp.cam.ac.uk
%

% returns eigenvalues of the laplacian finite difference matrix G for
% Neumann boundary conditions in matrix form, so for
%
%  Laplace(image) = reshape(G * image(:),size(image))
%
%  we have (for MATLAB's built-in dct functions or the equivalent dctt)
%
%    ev = laplace_dct_factors_neumann(size(image))
%    Laplace(image) = idct2(ev .* dct2(image))
%
%  To calculate the backward step (Laplace(image) + c)^-1 * b, use
%
%    x = idctN(dctN(image) ./ (ev + c)));
%
%  NOTE: we have
%    dct2(x)  = dct (dct (x )')',
%    idct2(x) = idct(idct(x')') ,
%  so for the backward step above:
%    x = idct(idct(((dct(dct(image)')') ./ (ev + c))')');
%
% "scaling" can be set (or set to [], or left out) to scale the gradient
% along the i-th dimension by scaling(i), i.e. A_i = A_i_base * scaling(i).
% Usually, set scaling = [1/h_x; 1/h_y; ...].

function ev = laplace_dct_factors_neumann(dims, scaling)
  if (nargin < 2 || isempty(scaling))
      scaling = ones(size(dims));
  end

  ev = zeros(prod(dims),1);
 
  for j = 1:numel(dims)
      kj = (0:(dims(j)-1))';
      evj = (1 - cos(3/2 * kj * pi/dims(j)) ./ cos(1/2 * kj * pi/dims(j))) .* (scaling(j)^2);

      ev = ev + kron(ones(prod(dims((j+1):end)),1),...
                kron(evj,ones(prod(dims(1:(j-1))),1)));
  end
  
  ev = reshape(ev,dims);

%{
  k1 = (0:(dims(1)-1))';
  ev1 = 1 - cos(3/2 * k1 * pi/dims(1)) ./ cos(1/2 * k1 * pi/dims(1));
  
  k2 = (0:(dims(2)-1))';
  ev2 = 1 - cos(3/2 * k2 * pi/dims(2)) ./ cos(1/2 * k2 * pi/dims(2));
  
  evx = reshape(kron(ones(dims(2),1),ev1) + kron(ev2,ones(dims(1),1)), dims);
  %  norm(ev(:)-evx(:),+inf)
%}  
end
