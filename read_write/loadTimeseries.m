%*******Alex Bhogal 1/29/2019
%function to load niftii images untouched. Output is the image data
%and info structure containing matrix and resolution info

function [image info] = loadTimeseries(pathname, filename)

warning('off')
global opts
if isfield(opts,'verbose'); else opts.verbose = 0; end
if ~ispc
    if filename(end) == char(10); filename(end) = []; end
end

image_path = [pathname,filename];
info = niftiinfo(image_path);
image = niftiread(image_path);


        if size(image,4) > 1
            if opts.verbose; disp('initializing opts.TR in opts structure - please check if correct otherwise modify opts.TR');
                disp('.');disp('.'); 
            end
            opts.TR = info.raw.pixdim(5);
            if opts.TR > 10; opts.TR = opts.TR/1000; end
            opts.dyn = size(image,4);
        end
        
        if isfield(opts,'info'); else
            if opts.verbose 
                disp('initializing timeseries header in opts.header.ts - used for saving 4D data');
                disp('.');disp('.');
                disp('initializing parameter map header in opts.header.map - used for saving 3D data');
                disp('.');disp('.');
                disp('initializing mask header in opts.header.map - used for saving masks');
                disp('.');disp('.'); 
            end
            opts.info.ts = info;
            opts.info.ts.Datatype = 'double';
            opts.info.ts.BitsPerPixel = 64
            opts.info.ts.raw.datatype = 64;
            opts.info.ts.raw.bitpix = 64;
            opts.info.map = info;
            opts.info.map.Datatype = 'double';
            opts.info.map.BitsPerPixel = 64
            opts.info.map.raw.datatype = 64;
            opts.info.map.raw.bitpix = 64;
            opts.info.map.raw.dim(1) = 3;
            opts.info.map.raw.dim(5) = 1;
            opts.info.map.raw.pixdim(5) = 0;
            opts.info.map.PixelDimensions = opts.info.ts.PixelDimensions(1:3);
            opts.info.map.ImageSize = opts.info.ts.ImageSize(1:3);
            %generate mask info
            opts.info.mask = opts.info.map;
            opts.info.mask.Datatype = 'double';
            opts.info.mask.BitsPerPixel = 64;
            opts.info.mask.raw.datatype = 64;
            opts.info.mask.raw.bitpix = 64;
            opts.voxelsize = info.raw.pixdim(2:4);
        end
 

image(isinf(image)) = 0;
image = double(image);
opts.niiwrite =1;
end


