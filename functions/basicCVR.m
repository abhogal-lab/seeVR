function [delta] = basicCVR(BOLD, WBmask, base_idx, stim_idx, savefile, opts)
%written by Alex Bhogal a.bhogal@umcutrecht.nl

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

