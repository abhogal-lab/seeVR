% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% sections of this code were contributed by Jeroen Siero, j.c.w.siero@umcutrecht.nl
% <smthData: smoothes input images using a gaussian kernel >
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
% smthData applies a 2D or 3D smoothing kernel to the input data.
%
% data: Timeseries data (i.e. 3D or 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.spatialdim, opts.voxelsize, opts.filtWidth, opts.FWHM
%
% sdata: this is smoothed version of the input data
function [sdata] = smthData(data, mask, opts)
mask = logical(mask);

warning('off')
global opts

if isfield(opts,'spatialdim'); else; opts.spatialdim = 2; end
if isfield(opts,'filtWidth'); else; opts.filtWidth = 5; end %originally 7
if isfield(opts,'FWHM'); else; opts.FWHM = 4; end
if isfield(opts,'verbose'); else; opts.verbose = 0; end %turn on/off select command output

sigma=opts.FWHM/2.355;
inplanevoxelsize=opts.voxelsize(1);
data_smooth=zeros(size(data));

filtSigma = sigma/inplanevoxelsize;

ndim = ndims(data);
%create a dilation mask
SE = strel('disk', 3);
dmask = imdilate(mask, SE);
delmask = double(dmask - mask); %cuts out center
ddata = zeros(size(data));

switch ndim
    case 3
        
        data(isnan(data)) = 0;
        ddata = imdilate(data,SE);% dilate dataset to pad
        
        ddata = ddata.*delmask;
        ddata = ddata + data.*mask;
        
        switch opts.spatialdim
            case 2
                imageFilter=fspecial('gaussian',opts.filtWidth,filtSigma);
                for s=1:size(data, 3)
                    if ~isempty(isnan(data))
                        
                        data_smooth(:,:,s)= nanconv(data(:,:,s),imageFilter, 'nanout');
                    else
                        data_smooth(:,:,s)=imgaussfilt(data(:,:,s),sigma/inplanevoxelsize, 'FilterDomain', 'spatial'); % here also fitler width of 7 is used internally
                    end
                end
            case 3
                if isnan(data)
                    
                    ddata(isnan(ddata)) = 0;
                    data_smooth=imgaussfilt3(ddata,sigma./opts.voxelsize, 'FilterDomain', 'spatial');
                else
                    data_smooth=imgaussfilt3(ddata,sigma./opts.voxelsize, 'FilterDomain', 'spatial');
                    
                end
        end
    case 4
        
        data(isnan(data)) = 0;
        
        for ii = 1:size(ddata,4)
            ddata(:,:,:,ii) = imdilate(data(:,:,:,ii),SE);% dilate dataset to pad
        end
        ddata = ddata.*delmask;
        ddata = ddata + data.*mask;
        
        switch opts.spatialdim
            case 2
                imageFilter=fspecial('gaussian',opts.filtWidth,filtSigma);
                
                for t=1:size(data,4)
                    for s=1:size(data, 3)
                        data_smooth(:,:,s,t)= nanconv(data(:,:,s,t),imageFilter, 'nanout');
                    end
                end
            case 3
                
                if isnan(data)
                    ddata(isnan(ddata)) = 0;
                    for t=1:size(data,4)
                        data_smooth(:,:,:,t)=imgaussfilt3(ddata(:,:,:,t),sigma./opts.voxelsize, 'FilterDomain', 'spatial');
                    end
                else
                    %dilate image to mitigate edge effects
                    for t=1:size(data,4)
                        data_smooth(:,:,:,t)=imgaussfilt3(ddata(:,:,:,t),sigma./opts.voxelsize, 'FilterDomain', 'spatial');
                    end
                end
        end
end
sdata = double(mask.*data_smooth);

disp('Finished smoothing Data')
end
