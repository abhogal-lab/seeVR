%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [idx, rdata] = chopTimeseries(data,mask,idx)

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

