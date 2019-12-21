function result = dimensions(imagedim, components)
    if (nargin < 2)
        components = size(imagedim, ndims(imagedim));
        imagedim = size(imagedim);
        imagedim = imagedim(1:end-1);
    end

    result = [];
    result.image = imagedim(:)';
    result.nimage = prod(imagedim);
    result.mimage = numel(imagedim);
    result.components = components;
    result.ntotal = result.nimage * result.components;
    result.total = [result.image components];
    result.mtotal = numel(result.total);
end
