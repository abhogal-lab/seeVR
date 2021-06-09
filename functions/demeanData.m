function [dmData] = demeanData(data, mask)
%written by Alex Bhogal a.bhogal@umcutrecht.nl
%function to generate mean timeseries from data based on input mask

disp('de-meaning timeseries data')
[x,y,z,dyn] = size(data);
[voxel_ts,coordinates] = grabTimeseries(data,mask);
m_ts = nanmean(voxel_ts,2);
voxel_ts = voxel_ts-m_ts;
dmData = zeros(size(data));
dmData = reshape(dmData,[x*y*z dyn]);
dmData(coordinates,:) = voxel_ts;
dmData = reshape(dmData,size(data));
end