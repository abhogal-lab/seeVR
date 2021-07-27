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

function [mmask] = remLV(data,mask,opts)
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

warning('off')
 global opts
           			
			[voxel_ts, coordinates] = grabTimeseries(data, mask);
            SD = nanstd(voxel_ts,0,2);
            MN = nanmean(voxel_ts,2);
            SDm = zeros(1,numel(mask)); tSNR = zeros(1,numel(mask));
            SDm(coordinates) = SD;
            SDm = reshape(SDm, size(mask));
            tSNR(coordinates) = MN./SD;
            tSNR = reshape(tSNR, size(mask));
            SDinv = 1/SDm;
            tNSR = 1/tSNR;
            
            %save image
            saveImageData(SDm,opts.headers.map,opts.resultsdir,'tSD.nii.gz',64);
            disp('saving tSD map')
            saveImageData(tSNR,opts.headers.map,opts.resultsdir,'tSNR.nii.gz',64);
            disp('saving tSNR map')
            saveImageData(SDinv,opts.headers.map,opts.resultsdir,'SDinv.nii.gz',64);
            disp('saving 1/tSD map')
            saveImageData(tNSR,opts.headers.map,opts.resultsdir,'tNSR.nii.gz',64);
            disp('saving tNSR map')
            
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
            
            disp('updating whole brain mask')
            mmask(isinf(tNSR)) = 0; mmask(tNSR > opts.LVthresh) = 0; %This step should remove vessels and garbage data
            mmask = double(mmask);
            %saves new brain mask excluding voxels
            saveImageData(mmask,opts.headers.mask,opts.resultsdir,['mWBmask_',num2str(opts.LVthresh),'.nii.gz'],64);

end

