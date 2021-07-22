%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [ndata] = normTimeseries(data,mask,idx)
%written by Alex Bhogal a.bhogal@umcutrecht.nl
%function to normalize timeseries to baseline using input mask
%index are some user-supplied baseline timepoints. If index is not empty, this
%function will allow manual selection of baseline timepoints
global opts;

p = cputime;
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
disp(['Dataset has been normalized to baseline period']);
end