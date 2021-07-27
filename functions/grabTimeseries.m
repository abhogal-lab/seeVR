%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [voxel_ts, coordinates] = grabTimeseries(data, mask)

switch ndims(data)
    case 4
[x,y,z,dyn] = size(data); 
mask = logical(mask(:));
coordinates = find(mask); 
%disp('grabbing timeseries (4D) at each voxel')

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

    case 3
[x,y,dyn] = size(data);
mask = logical(mask(:));
coordinates = find(mask); 
%disp('grabbing timeseries (3D) at each voxel')

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

end
end