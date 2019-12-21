
function result = projector_uniform_unitsimplex(dims)
    result = @(u)(reshape(...
        project_unitsimplex(reshape(u,[dims.nimage dims.components])),...
        [dims.ntotal 1]));
end