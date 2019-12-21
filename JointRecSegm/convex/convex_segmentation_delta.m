

function [f, fd, u, ud]=convex_segmentation_delta(image, lambda, classes, delta)

    dbglevel(5);
    rand('seed',192497); % make noise predictable

 
    f = zeros(size(image));
    for i = 1:size(classes,1)
        f(:,:,i) =  delta * (image - classes(i)).^2; 
    end

    % Set up the problem for the solver
    %magedim, components
    dims = dimensions([size(f,1) size(f,2)],size(classes,1));
    
    
    p = 2; % which norm to use for the regularizer (always use 2 as it is isotropic)
    
    %tv = regularizer_uniform_tv_fast(dims, lambda, p, 'forward','neumann'); % this requires to compile a MEX file
    tv = regularizer_uniform_tv(dims, lambda, p, 'forward', 'neumann');
    unitsimplex = constraints_uniform_unitsimplex(dims);
    dataterm = dataterm_linear(f(:));

%     info = @(u,i,d,details)print_info(u,i,d,details, 1,dataterm,unitsimplex,tv);
%     ppd = struct(...
%         'term_maxiter',10000,...
%         'term_optimality_relative',1e-6,...        
%         'callback_details',info);
%     
    % solve the problem
    
    [u,~] = fpd(unitsimplex.center, 0*tv.dual_constraints.center, dataterm, unitsimplex, tv, ppd);
    
    % discretise u
    
    ud = discretize_firstmax(u, dims);
    
    % reshape u and ud 
    
    u = reshape(u, size(f));
    ud = reshape(ud, size(f));
    
    fd = discretize_firstmax(-f, dims);
    fd = reshape(fd, size(f));

    
end
