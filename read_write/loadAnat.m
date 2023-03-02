%*******Alex Bhogal 1/29/2019
%function to load niftii images untouched. Output is the image data
%and info structure containing matrix and resolution info

function [image info] = loadAnat(pathname, filename)

warning('off')
global opts;
if isfield(opts,'verbose'); else opts.verbose = 0; end
if ~ispc
    if filename(end) == char(10); filename(end) = []; end
end

image_path = fullfile(pathname,filename);
info = niftiinfo(image_path);
image = niftiread(image_path);
opts.anatpath = pathname;
opts.anatfile = filename;
opts.voxelsize_anat = info.PixelDimensions;

%generate anat info
opts.info.anat = info;
opts.info.anat.Datatype = 'double';
opts.info.anat.BitsPerPixel = 64;
opts.info.anat.raw.datatype = 64;
opts.info.anat.raw.bitpix = 64;
opts.info.anat.MultiplicativeScaling = 1;
image(isinf(image)) = 0;
image = double(image);
opts.niiwrite = 1;

end


