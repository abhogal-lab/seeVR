function [voxel_ts, coordinates] = grabTimeseries(data, mask)
%written by Alex Bhogal a.bhogal@umcutrecht.nl
switch ndims(data)
    case 4
[x,y,z,dyn] = size(data); 
mask = logical(mask(:));
coordinates = find(mask); 
disp('grabbing timeseries (4D) at each voxel')
tic
clean_voxel_ts=zeros([length(coordinates) dyn]);
data = reshape(data,[x*y*z dyn]);
clean_voxel_ts = data(coordinates,:);
clean_voxel_ts(clean_voxel_ts == 0) = NaN;
voxel_ts = clean_voxel_ts;
%remove NaN rows
voxel_ts(any(isnan(clean_voxel_ts), 2), :) = [];
coordinates(any(isnan(clean_voxel_ts), 2), :) = [];
%remove Inf rows
coordinates(any(isinf(voxel_ts), 2), :) = [];
voxel_ts(any(isinf(voxel_ts), 2), :) = [];
toc
    case 3
[x,y,dyn] = size(data);
mask = logical(mask(:));
coordinates = find(mask); 
disp('grabbing timeseries (3D) at each voxel')
tic
clean_voxel_ts=zeros([length(coordinates) dyn]);
data = reshape(data,[x*y dyn]);
clean_voxel_ts = data(coordinates,:);
clean_voxel_ts(clean_voxel_ts == 0) = NaN;
voxel_ts = clean_voxel_ts;
%remove NaN rows
voxel_ts(any(isnan(clean_voxel_ts), 2), :) = [];
coordinates(any(isnan(clean_voxel_ts), 2), :) = [];
%remove Inf rows
coordinates(any(isinf(voxel_ts), 2), :) = [];
voxel_ts(any(isinf(voxel_ts), 2), :) = [];
toc
end
end