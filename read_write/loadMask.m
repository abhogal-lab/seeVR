%*******Alex Bhogal 1/29/2019
%function to load niftii images untouched. Output is the image data
%and info structure containing matrix and resolution info

function [image info] = loadMask(pathname, filename)

warning('off')
global opts;
if isfield(opts,'verbose'); else opts.verbose = 0; end
if ~ispc
    if filename(end) == char(10); filename(end) = []; end
end

image_path = [pathname,filename];
info = niftiinfo(image_path);
image = niftiread(image_path);

            %generate mask info
            opts.info.mask = info;
            opts.info.mask.Datatype = 'double';
            opts.info.mask.BitsPerPixel = 64;
            opts.info.mask.raw.datatype = 64;
            opts.info.mask.raw.bitpix = 64;
            opts.info.mask.MultiplicativeScaling = 1;
        image(isinf(image)) = 0;
image = double(image);
opts.niiwrite = 1; 

end


