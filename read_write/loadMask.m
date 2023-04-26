function [image info] = loadMask(pathname, filename)

warning('off')
global opts;
if isfield(opts,'verbose'); else opts.verbose = 0; end
if ~ispc
    if filename(end) == char(10); filename(end) = []; end
end

image_path = fullfile(pathname,filename);
info = niftiinfo(image_path);
image = niftiread(image_path);

%generate mask info
opts.info.mask = info;
opts.info.mask.MultiplicativeScaling = 1;
opts.niiwrite = 1;
opts.maskDatatype = info.Datatype;

end


