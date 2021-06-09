
function [ndata] = normTimeseries(data,mask,idx)
%written by Alex Bhogal a.bhogal@umcutrecht.nl
%function to normalize timeseries to baseline using input mask
%data is the timeseries data to be normalized
%mask is any corresponding mask to isolate voxels of interest
%index are some user-supplied baseline timepoints. If index is not empty, this
%function will allow manual selection of baseline timepoints
global opts;

p = cputime;
switch nargin
    case 2
TS = meanTimeseries(data,mask);
figure(100); plot(TS); title('Isolate first baseline period for normalization: 2 clicks'); 
[idx,~] = ginput(2); idx = round(idx); %select points 
close(100);
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
disp(['Dataset was normalized to baseline period in ',int2str(cputime-p),' seconds']);
end