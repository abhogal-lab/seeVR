function [img_ts] = meanTimeseries(data, mask)
%written by Alex Bhogal a.bhogal@umcutrecht.nl
%function to generate mean timeseries from data based on input mask

mask(isnan(mask)) = 0;
mask = logical(mask(:));
coordinates = find(mask); 

disp('Generate mean timeseries')

[voxel_ts,~] = grabTimeseries(data,mask);

img_ts = nanmean(voxel_ts,1);   
end