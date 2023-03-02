% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <grabTimeseres: extract voxel timeseries data and associated coordinates >
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% *************************************************************************
% takes 3D or 4D input data and extracts voxel timeseries data defined by
% the input mask
%
% data: input timeseries data (i.e. 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% voxel_ts: a tall array with each voxel timeseries defined by the input
% mask
%
% coordinates: a tall array of spatial coordinates associated with each
% voxel timeseries
function [voxel_ts, coordinates] = grabTimeseries(data, mask)
tf = class(data); 

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
voxel_ts = cast(voxel_ts,tf); coordinates = cast(coordinates,tf);

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
voxel_ts = cast(voxel_ts,tf); coordinates = cast(coordinates,tf);
end
end