%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [img_ts] = meanTimeseries(data, mask)
%function to generate mean timeseries from data based on input mask

mask(isnan(mask)) = 0;
mask = logical(mask(:));
coordinates = find(mask); 

disp('Generate mean timeseries')

[voxel_ts,~] = grabTimeseries(data,mask);

img_ts = nanmean(voxel_ts,1);   
end