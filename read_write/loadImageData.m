%*******Alex Bhogal 1/29/2019
%function to load niftii images untouched. Output is the image data
%and info structure containing matrix and resolution info

function [image info header] = loadImageData(pathname, filename)

warning('off')
global opts
if isfield(opts,'verbose'); else opts.verbose = 0; end;
if ~ispc
    if filename(end) == char(10); filename(end) = []; end
end

image_path = fullfile(pathname,filename);
data = load_untouch_nii(image_path);
image = (data.img);
info.datatype = data.hdr.dime.datatype;
info.xdim = data.hdr.dime.dim(2);
info.ydim = data.hdr.dime.dim(3);
info.zdim = data.hdr.dime.dim(4);
info.resx = data.hdr.dime.pixdim(2);
info.resy = data.hdr.dime.pixdim(3);
info.resz = data.hdr.dime.pixdim(4);

info.dim = [info.resx info.resy info.resz];

header = load_untouch_header_only(image_path);

switch ndims(image)
    case 4
        if header.dime.dim(5) > 1
            if opts.verbose; disp('initializing opts.TR in opts structure - please check if correct otherwise modify opts.TR');
                disp('.');disp('.'); 
            end
            opts.TR = header.dime.pixdim(5);
            if opts.TR > 10; opts.TR = opts.TR/1000; end
            opts.dyn = size(image,4);
            opts.xdata = (opts.TR:opts.TR:opts.TR*opts.dyn);
        end
        if isfield(opts,'headers'); else
            if opts.verbose; disp('initializing timeseries header in opts.header.ts - used for saving 4D data');
                disp('.');disp('.');
                disp('initializing parameter map header in opts.header.map - used for saving 3D data');
                disp('.');disp('.');
                disp('initializing mask header in opts.header.map - used for saving masks');
                disp('.');disp('.'); 
            end
            opts.headers.ts = header;
            opts.headers.map =  header;
            opts.headers.map.dime.dim(5) = 1;
            opts.headers.mask = opts.headers.map;
            opts.headers.mask.dime.datatype = 4;
            opts.voxelsize = header.dime.pixdim(2:4);
        end
    case 3
        if header.dime.datatype == 4
            if opts.verbose; disp('Based on the datatype this could be a mask')
                disp('.');disp('.'); 
            end
            opts.headers.mask = header;
            opts.headers.mask.dime.dim(5) = 1;
            opts.headers.mask.dime.datatype = 4;
            opts.voxelsize = header.dime.pixdim(2:4);
        else
            if opts.verbose; disp('initializing parameter map header in opts.header.map - used for 3D data');
                disp('.');disp('.'); 
            end
            opts.headers.map =  header;
            opts.headers.map.dime.dim(5) = 1;
            opts.voxelsize = header.dime.pixdim(2:4);
            if isfield(opts.headers,'mask')
                if opts.verbose; disp('Mask header information already exists in opts.headers.mask'); end
            else
                if opts.verbose; disp('Initializing mask header in opts.header.mask - used for saving masks');
                    disp('.');disp('.'); 
                end
                opts.headers.mask = header;
                opts.headers.mask.dime.dim(5) = 1;
                opts.headers.mask.dime.datatype = 4;
                opts.voxelsize = header.dime.pixdim(2:4);
            end
        end
end

image(isinf(image)) = 0;
image = double(image);
end


