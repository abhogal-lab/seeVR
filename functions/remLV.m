%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [mmask] = remLV(data,mask,opts)
%Written by Alex Bhogal, a.bhogal@umcutrecht.nl
%This function uses the inverse tSNR (tNSR) to isolate bright signals associated
%with large draining veins. These voxels are removed from the original
%whole brain mask and returned for further processing: 1/tSD, tSD, tSNR, tNSR are
%saved
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

