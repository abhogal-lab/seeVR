%%%ALEX BHOGAL
%function that displays timeseries data as a quick check
%data = 4D timeseries
%slice = which slice of the 3D volume to show
%scale = range of the data to be shown

function [] = view_ts(data,slice,scale)

figure;
for ii=1:size(data,4)
    imagesc(imrotate(data(:,:,slice,ii),90), scale)
    colormap(jet)
    pause(0.01)
end

end