%*******Alex Bhogal 1/29/2019
%function to load niftii images untouched. Output is the image data
%and info structure containing matrix and resolution info

function [image info] = loadAnat(pathname, filename)

warning('off')
global opts;

if ~ispc
    if filename(end) == char(10); filename(end) = []; end
end

image_path = fullfile(pathname,filename);
info = niftiinfo(image_path);
image = niftiread(image_path);
opts.anatpath = pathname;
opts.anatfile = filename;
opts.voxelsize_anat = info.PixelDimensions;
opts.anatDatatype = info.Datatype;
%generate anat info
opts.info.anat = info;
image(isinf(image)) = 0;
opts.niiwrite = 1;

end


