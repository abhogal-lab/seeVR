% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <remLV: removes unwanted signals from large vessels or CSF >
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
% remLV function uses the inverse tSNR (tNSR) to isolate bright signals associated
% with large draining veins. These voxels are removed from the original
% mask and returned for further processing: 1/tSD, tSD, tSNR, tNSR are
% saved
%
% data: input timeseries data (i.e. 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.LVpercentile OR opts.LVthresh, opts.headers.map, opts.resultsdir
%
% mmask: this is the modified mask where voxels above the specified
% threshold have been removed
function [mmask] = remLV(data,mask,opts)

warning('off')
global opts

if isfield(opts,'verbose'); else; opts.verbose = 0; end %turn on/off select command output
if isfield(opts,'niiwrite'); else; opts.niiwrite = 0; end


if isfield(opts,'resultsdir'); else; opts.resultsdir = fullfile(pwd); end

[x, y, z, ~] = size(data);
[voxel_ts, coordinates] = grabTimeseries(data, mask);
SD = nanstd(voxel_ts,0,2);
MN = nanmean(voxel_ts,2);
SDm = zeros(1,x*y*z); tSNR = zeros(1,x*y*z);
SDm(coordinates) = SD;
SDm = reshape(SDm, size(mask));
tSNR(coordinates) = MN./SD;
tSNR = reshape(tSNR, size(mask));
SDinv = 1/SDm;
tNSR = 1/tSNR;

%save image

if opts.niiwrite
    cd(opts.resultsdir);
    niftiwrite(cast(SDm,opts.mapDatatype),'tSD',opts.info.map);
    niftiwrite(cast(tSNR,opts.mapDatatype),'tSNR',opts.info.map);
    niftiwrite(cast(SDinv,opts.mapDatatype),'SDinv',opts.info.map);
    niftiwrite(cast(tNSR,opts.mapDatatype),'tNSR',opts.info.map);
else
    saveImageData(SDm,opts.headers.map,opts.resultsdir,'tSD.nii.gz',64);
    saveImageData(tSNR,opts.headers.map,opts.resultsdir,'tSNR.nii.gz',64);
    saveImageData(SDinv,opts.headers.map,opts.resultsdir,'SDinv.nii.gz',64);
    saveImageData(tNSR,opts.headers.map,opts.resultsdir,'tNSR.nii.gz',64);
end

%update WBmask to exclude unreliable signals
mmask = mask;
tNSR(isinf(tNSR)) = 0;
tNSR(isnan(tNSR)) = 0;

if isfield(opts,'LVpercentile')
    opts.LVthresh = prctile(tNSR(:),opts.LVpercentile);
elseif isfield(opts,'LVthresh')
else
    opts.LVthresh = prctile(tNSR(:),98);
end

mmask(isinf(tNSR)) = 0; mmask(tNSR > opts.LVthresh) = 0; %This step should remove vessels and garbage data
%saves new brain mask excluding voxels
if opts.niiwrite
    name = fullfile(opts.resultsdir,['mWBmask_',num2str(opts.LVthresh)]);
    name(name == '.') = '_';
    niftiwrite(cast(mmask,opts.maskDatatype),name,opts.info.mask);
else
    saveImageData(mmask,opts.headers.mask,opts.resultsdir,['mWBmask_',num2str(opts.LVthresh),'.nii.gz'],64);
end

    disp('saving tSD map');
    disp('saving tSNR map');
    disp('saving 1/tSD map')
    disp('saving tNSR map')
    disp('updating whole brain mask')
    
    % Plot
    figure;
    z = round(size(data,3)/2);
    subplot(1,3,1)
    imagesc(imrotate(mean(data(:,:,z,:),4),90));
    title('mean image', 'Color', 'w'); colormap(gray); cleanPlot;
    subplot(1,3,2)
    imagesc(imrotate(mask(:,:,z),90));
    title('original mask', 'Color', 'w'); cleanPlot;
    subplot(1,3,3)
    imagesc(imrotate(mmask(:,:,z),90));
    title('modified mask', 'Color', 'w'); cleanPlot;
    set(gcf, 'Units', 'pixels', 'Position', [200, 500, 600, 160]);

end

