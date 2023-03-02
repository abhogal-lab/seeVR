% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <normTimeseres: calculates %change from baseline signal >
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
% This function will normalize timeseries data (i.e. calculate %change) to
% baseline using input mask index are some user-supplied baseline
% timepoints. If index is not empty, this function will allow manual 
% selection of baseline timepoints
%
% data: input timeseries data (i.e. 3D (slice) or 4D (volume) MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% idx: this can be either a 2 element vector consisting of start and end
% points. If this vector is not supplied, the user will be promted to
% select the beginning and ending of the desired epoch manually. 
%
% ndata: the baseline normalized version of the input data
%
% opts.norm_idx: this parameter is changed to reflect the selected baseline
% start and end indices.
function [ndata] = normTimeseries(data,mask,idx)

warning('off')
global opts;
tf = class(data);
switch nargin
    case 2
TS = meanTimeseries(data,mask);
figure; plot(TS); title('Select start and end of baseline: 2 clicks'); 
[idx,~] = ginput(2); idx = round(idx); %select points 
close(gcf);
    case 3
end

[voxel_ts coordinates] = grabTimeseries(data, mask);
bvoxel_ts = mean(voxel_ts(:,idx(1):idx(2)),2);

%normalize signal
switch ndims(data)
    case 3
[x,y,dyn] = size(data);     
voxel_ts = 100*((voxel_ts - repmat(bvoxel_ts,1,dyn))./bvoxel_ts);
%re-fill normalized BOLD timeseries
ndata = nan([x*y dyn]);
ndata(coordinates,:) = voxel_ts;
ndata = reshape(ndata, [x y dyn]);
case 4
[x,y,z,dyn] = size(data);     
voxel_ts = 100*((voxel_ts - repmat(bvoxel_ts,1,dyn))./bvoxel_ts);
%re-fill normalized BOLD timeseries
ndata = nan([x*y*z dyn]); 
ndata(coordinates,:) = voxel_ts;
ndata = reshape(ndata, [x y z dyn]);
end
opts.norm_idx = idx;
%disp(['Dataset has been normalized to baseline period']);
ndata = cast(ndata,tf);
end