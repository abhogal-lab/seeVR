function [idx, rdata] = chopTimeseries(data,mask,idx)
%written by Alex Bhogal a.bhogal@umcutrecht.nl
%function to generate mean timeseries from data based on input mask
%data is the timeseries data to be shortened
%mask is any corresponding mask to isolate voxels of interest
%index are some user-supplied timepoints. If index is not empty, this
%function will allow manual selection of timepoints
global opts
switch nargin
    case 2
        figure;  plot(meanTimeseries(data,mask));
        title('Select start and end point of epoch');
        [idx,~] = ginput(2); idx = round(idx);
        close;
        switch ndims(data)
            case 4
                rdata = data(:,:,:,idx(1):idx(2));
            case 3
                rdata = data(:,:,idx(1):idx(2));
        end
        
    case 3
        switch ndims(data)
            case 4
                rdata = data(:,:,:,idx(1):idx(2));
            case 3
                rdata = data(:,:,idx(1):idx(2));
        end
        disp('Using user-supplied indices to select epoch')
end
disp('Updated header file to reflect new timeseries lenght')
opts.headers.ts.dime.dim(2:5) = size(rdata);
opts.dyn = []; [opts.xdim,opts.ydim,opts.zdim,opts.dyn] = size(rdata);
opts.xdata = [opts.TR:opts.TR:opts.TR*opts.dyn];
disp('Updated opts.xdata and opts.dyn to reflect new timeseries length');
end

