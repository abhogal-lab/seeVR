function [img info] = loadImage(pathname, filename)

warning('off')
global opts;
if isfield(opts,'verbose'); else opts.verbose = 0; end
if ~ispc
    if filename(end) == char(10); filename(end) = []; end
end

image_path = fullfile(pathname,filename);
info = niftiinfo(image_path);
img = niftiread(image_path);
opts.imagepath = pathname;
opts.imagefile = filename;
opts.voxelsize_image = info.PixelDimensions;
opts.imageDatatype = info.Datatype;

%generate anat info
opts.info.image = info;
img(isinf(img)) = 0;
opts.niiwrite = 1;

end


