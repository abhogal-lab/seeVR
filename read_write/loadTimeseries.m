function [image info] = loadTimeseries(pathname, filename)

warning('off')
global opts;
if isfield(opts,'verbose'); else opts.verbose = 0; end
if ~ispc
    if filename(end) == char(10); filename(end) = []; end
end

image_path = fullfile(pathname,filename);
info = niftiinfo(image_path);
image = niftiread(image_path);
opts.funcpath = pathname;
opts.funcfile = filename;
opts.voxelsize = info.PixelDimensions(1:3);

if size(image,4) > 1
    if opts.verbose; disp('initializing opts.TR in opts structure - please check if correct otherwise modify opts.TR');
        disp('.');disp('.');
    end
    opts.TR = info.raw.pixdim(5);
    if opts.TR > 10; opts.TR = opts.TR/1000; end
    opts.dyn = size(image,4);
end


opts.info.ts = info;
opts.tsDatatype = info.Datatype;

opts.info.map = info;
opts.info.map.raw.dim(1) = 3;
opts.info.map.raw.dim(5) = 1;
opts.info.map.raw.pixdim(5) = 0;
opts.info.map.PixelDimensions = opts.info.ts.PixelDimensions(1:3);
opts.info.map.ImageSize = opts.info.ts.ImageSize(1:3);
opts.mapDatatype = info.Datatype;

image(isinf(image)) = 0;

opts.niiwrite = 1;
end


