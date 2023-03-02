% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <meanTimeseres:  generates a mean timeseries from data based on input mask>
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
% function to generate mean timeseries from data based on input mask
%
% data: input timeseries data (i.e. 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
function [img_ts] = meanTimeseries(data, mask)

mask(isnan(mask)) = 0;
mask = logical(mask(:));
coordinates = find(mask); 
[voxel_ts,~] = grabTimeseries(data,mask);
img_ts = nanmean(voxel_ts,1);   

end