
%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [delta] = basicCVR(BOLD, WBmask, base_idx, stim_idx, savefile, opts)
if isfield(opts,'smoothmap'); else; opts.smoothmap = 1; end;

disp('Normalizing to baseline')
BOLD = normTimeseries(BOLD,WBmask,base_idx);
sBOLD = BOLD; sBOLD(isnan(sBOLD)) = 0;

sBOLD( sBOLD == 0) = NaN;
BOLD = sBOLD; clear sBOLD;

base = nanmean(BOLD(:,:,:,base_idx(1):base_idx(2)),4);
stim = nanmean(BOLD(:,:,:,stim_idx(1):stim_idx(2)),4);
delta = stim - base; delta(~WBmask) = NaN;

%spatially smooth map (only done if timeseries is not smoothed already)
if opts.smoothmap
delta = smthData(delta , WBmask, opts);
end
delta(delta > 15) = 0; delta(delta < -15) = 0; %remove large BOLD
delta(isnan(delta)) = 0;
disp('Saving CVR map')
saveImageData(delta,opts.headers.map,opts.resultsdir,savefile, 64);
end

