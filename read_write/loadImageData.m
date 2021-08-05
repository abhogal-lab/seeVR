%*******Alex Bhogal 1/29/2019
%function to load niftii images untouched. Output is the image data
%and info structure containing matrix and resolution info

function [image info header] = loadImageData(pathname, filename)

if ~ispc
  if filename(end) == char(10); filename(end) = []; end  
end
image_path = [pathname,filename];
data = load_untouch_nii(image_path);
image = double(data.img);
info.datatype = data.hdr.dime.datatype;
info.xdim = data.hdr.dime.dim(2);
info.xdim = data.hdr.dime.dim(3);
info.zdim = data.hdr.dime.dim(4);
if size(data.hdr.dime.dim,2)>4
info.dyn = data.hdr.dime.dim(5);
end

info.resx = data.hdr.dime.pixdim(2);
info.resy = data.hdr.dime.pixdim(3);
info.resz = data.hdr.dime.pixdim(4);

info.dim = [info.resx info.resy info.resz];

header = load_untouch_header_only(image_path);



end