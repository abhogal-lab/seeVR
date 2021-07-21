%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [] = view_ts(data,slice,scale)
%function that displays timeseries data as a quick check
%data = 4D timeseries
%slice = which slice of the 3D volume to show
%scale = range of the data to be shown
figure;
for ii=1:size(data,4)
    imagesc(imrotate(data(:,:,slice,ii),90), scale)
    colormap(jet)
    pause(0.01)
end

end