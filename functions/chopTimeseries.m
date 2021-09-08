% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <chopTimeseres: shortens timeseries data to within a specified epoch >
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
% chopTimeseries shortens the input data in the temporal dimension as
% specified by either user input or supplied indices.
%
% data: input timeseries data (i.e. 3D (slice) or 4D (volume) MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% idx: this can be either a 2 element vector consisting of start and end
% points. If this vector is not supplied, the user will be promted to
% select the beginning and ending of the desired epoch manually. For manual
% selection, the indices (idx) are returned
%
% rdata: this is the shortened version of the input data
%
% opts.xdata and opts.dyn are updated to reflect the new shortened
function [idx, rdata] = chopTimeseries(data,mask,idx)

% timeseries
warning('off');
global opts;

if isfield(opts,'save_rdata'); else; opts.save_rdata = 0; end %saves the shortened timeseries
if isfield(opts,'verbose'); else; opts.verbose = 0; end %turn on/off select command output
if isfield(opts,'dyn'); else; opts.dyn = size(data,ndims(data)); end %setup option if its not there
if isfield(opts,'xdata'); else; opts.xdata = [opts.TR:opts.TR:opts.TR*opts.dyn]; end %setup option if its not there yet


switch nargin
    case 2
        figure;  plot(meanTimeseries(data,mask));
        title('Select start and end point of epoch');
        [idx,~] = ginput(2); idx = round(idx);
        close;
        
        switch ndims(data)
            case 4
                if idx(2) > size(data,4)
                    idx(2) = size(data,4);
                end
                rdata = data(:,:,:,idx(1):idx(2));
            case 3
                if idx(2) > size(data,3)
                    idx(2) = size(data,3);
                end
                rdata = data(:,:,idx(1):idx(2));
        end
        
    case 3
        switch ndims(data)
            case 4
                if idx(2) > size(data,4)
                    idx(2) = size(data,4);
                end
                rdata = data(:,:,:,idx(1):idx(2));
            case 3
                 if idx(2) > size(data,3)
                    idx(2) = size(data,3);
                end
                rdata = data(:,:,idx(1):idx(2));
        end
        if opts.verbose; disp('Using user-supplied indices to select epoch'); end
end
if opts.verbose
    disp('Updated header file to reflect new timeseries length'); 
    disp('Updated opts.xdata and opts.dyn to reflect new timeseries length'); 
end

if isfield(opts,'headers.ts'); opts.headers.ts.dime.dim(2:5) = size(rdata); end
if isfield(opts,'dyn'); opts.dyn = []; [opts.xdim,opts.ydim,opts.zdim,opts.dyn] = size(rdata); end
if isfield(opts,'xdata'); opts.xdata = [opts.TR:opts.TR:opts.TR*opts.dyn]; end

%save shortened timeseries (this needs to be updated for native saving)
if opts.save_rdata
    if opts.niiwrite
    niftiwrite(rdata,[opts.glmlagdir, 'rBOLD',],'Compressed',1); %need to add info for compatibility
    else
    saveImageData(rdata, opts.headers.ts, opts.resultsdir, 'rBOLD.nii.gz', 64);
    end
end
end

