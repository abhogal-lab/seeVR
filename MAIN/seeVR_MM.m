%%written by Alex Bhogal (a.bhogal@umcutrecht.nl). Do not share without
%%explicit permission from A. Bhogal

clear all
close all

%load acz MOYAMOYA BOLD

%for subject = [3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25]
subject = [9] %select individual subject. Alternatively run for loop to process all.
if subject < 10
    subj = ['MM0',int2str(subject)]
else
    subj = ['MM',int2str(subject)]
end

dir = ['C:\Users\abhogal\Documents\DATA\MRI_data\MOYAMOYA\Subjects\',subj,'\'];
cd(dir)

%filename = ls('*_masked_mcf.nii*')
filename = ls('BOLD_applytopup.nii.gz*')
[sourceBOLD,INFO,BOLDheader] = loadImageData(dir, filename);

%load GM mask
maskname = ls('*seg_1*')
[WMmask,INFOmask,HEADERmask] = loadImageData(dir, maskname);
maskname = ls('*seg_0*')
[GMmask,~,~] = loadImageData(dir, maskname);
maskname = ls('*seg_2*')
[CSFmask,~,~] = loadImageData(dir, maskname);
maskname = ls('*mean_brain_mask*')
[WBmask,~,~] = loadImageData(dir, maskname);

%load sinus or cerebellum masks that you make
maskname = ls('*cerebellum*')
[CERmask,~,~] = loadImageData(dir, maskname);

%% %setup some basic data options
opts.TR = BOLDheader.dime.pixdim(5); if opts.TR > 10; opts.TR = opts.TR/1000; end;
opts.voxelsize = BOLDheader.dime.pixdim(2:4);
[opts.xdim, opts.ydim, opts.zdim, opts.dyn] = size(sourceBOLD);
xdata = [opts.TR:opts.TR:opts.TR*opts.dyn];

%%%setup some basic folder options
opts.savedir = dir;
opts.headers = struct();
opts.headers.ts = BOLDheader;
opts.headers.map = HEADERmask;
cd(opts.savedir)
opts.figdir = [dir,'figures\']; mkdir(opts.figdir)
opts.resultsdir = [dir,'RESULTS\'];  mkdir(opts.resultsdir)

%%%setup some preprocessing options
opts.regr_mp = 0;
opts.motioncorr = 0.7; %threshold to include motion regressors; for MM data these are often highly correlated with ACZ effect
opts.legOrder = [0]; %order of Legendre polynomials to include in regression (typically up to 4th order however can cause problems)
opts.plot = 0; %turn on plotting in functions that include this option
opts.LVthresh = 0.06; %threshold for removing vessel signals by remLV function (trial&error depending on data)
opts.smoothTS = 1; %spatially smooth BOLD timeseries
opts.smoothmap = 0; %certain functions will smoot 2D maps if this is on. For now, timeseries data is smoothed already so this is not needed.
%spatially smooth data
opts.spatialdim = 3;% 2 for 2D smoothing 3 for 3D. Edge effects for 3D are handled with dilation algorithm
opts.FWHM = 5; %ASL uses FWHM between 5-8mm; can consider to increase for comparison
opts.denoise = 1; %perform wavelet denoising if toolbox is available (otherwise temporal smoothing)

%% remove large vessel contributions, update WBmask, save tSNR/SD etc
[mWBmask] = remLV(sourceBOLD,WBmask,opts);
%temporally denoise data
%% denoise and smooth data

if opts.denoise; [BOLD] = denoiseData(sourceBOLD,WBmask,opts); end
if opts.smoothTS; BOLD = smthData( BOLD, WBmask, opts); end
%BOLD = sourceBOLD
%generate basic CVR map
stim_idx = [size(BOLD,4)-35 size(BOLD,4)];
base_idx = [5 20]; %indeces used to normalized BOLD data

delta = basicCVR(BOLD, WBmask, base_idx, stim_idx, ['delta_s_',subj,'_base.nii.gz'],opts);

n_mask = delta; n_mask(n_mask >-0.01) = 0; n_mask(n_mask <0) = 1;
p_mask = delta; p_mask(p_mask >0.01) = 1; p_mask(p_mask <0.01) = 0;

nTS = make_mean_ts(BOLD, n_mask); %mask of negative response regions
pTS = make_mean_ts(BOLD, p_mask); %mask of positive response regions

%% process MM data

%load motion parameters
mpfilename = ls('*mcf.par')
nuissance = load(mpfilename);
autoCorrP = abs(corr(pTS',nuissance)); autoCorrN = abs(corr(nTS',nuissance));
autoCorrP(autoCorrP < opts.motioncorr) = 0; autoCorrP(autoCorrP > 0) = 1; autoCorrN(autoCorrN < opts.motioncorr) = 0; autoCorrN(autoCorrN > 0) = 1;
% %remove highly correlated nuissance regressors to preserve signal response
autoCorr = autoCorrP + autoCorrP; autoCorr(autoCorr > 0 ) = 1;
index = ([1:1:size(nuissance,2)]).*autoCorr; index(index == 0) = [];
nuissance(:,index) = [];
dtnuissance = gradient(nuissance);
sqnuissance = nuissance.*nuissance;

%initialize Legendre regressors
if opts.legOrder == 0; L = []; else
    L = []; for ii=opts.legOrder; L(ii,:) = rescale(LegendreN(ii,xdata)); end
    L = L';
end

M =[nuissance dtnuissance sqnuissance]; %M = detrend(M); %can try detrending motion params for MM data
fd = [L M]; %combine nuisance regressors (i.e. Legendre, motion params) 

%fit a SLM model to positive series
%positive regions (from masks generated above)
slm = slmengine(double(xdata),double(pTS),'plot','on','knots',6,'increasing','on', ...
    'leftslope',0,'rightslope',0);  saveas(gcf,[opts.figdir,'positive_fit_',subj,'.fig']);
fpT = slmeval(xdata,slm)

%establise nuissance regressors
nuisPTS = pTS - fpT;
%negative regions (from masks generated above)
slm = slmengine(double(xdata),double(nTS),'plot','on','knots',6,'increasing','off', ...
    'leftslope',0,'rightslope',0);  saveas(gcf,[opts.figdir,'negative_fit_',subj,'.fig']);
fnT = slmeval(xdata,slm)

%establish nuissance regressors
nuisNTS = nTS - fnT;

%% Generate HRF
%positive response
opts.rratio = 6; %seems to affect the spread (use large value to limit)
opts.onset = 1:10*opts.TR:80*opts.TR;
opts.disp = 0.4;
opts.under = 1;
[HRF,onHRFprobe] = convHRF(fpT,opts); %Generate a double gamma HRF focus on onset
%negative response
[nHRF,nonHRFprobe] = convHRF(fnT,opts); %Generate a double gamma HRF focus on dispersion

cshift = -round(20/opts.TR);
for ii=1:size(onHRFprobe,1)
    onHRFprobe(ii,:) = circshift(onHRFprobe(ii,:),cshift); onHRFprobe(ii,end+cshift:end) = onHRFprobe(ii,end+cshift);
    nonHRFprobe(ii,:) = circshift(nonHRFprobe(ii,:),cshift); nonHRFprobe(ii,end+cshift:end) = nonHRFprobe(ii,end+cshift);
end
figure(30); subplot(1,2,1); plot([onHRFprobe' nonHRFprobe']); title('Shifted onset, constant dispersion')

%positive response
opts.onset = 1;
opts.disp = 1:10*opts.TR:80*opts.TR;
opts.under = 1;
[HRF,dspHRFprobe] = convHRF(fpT,opts); %Generate a double gamma HRF focus on onset
%negative response
[nHRF,ndspHRFprobe] = convHRF(fnT,opts); %Generate a double gamma HRF focus on dispersion
figure(30); subplot(1,2,2); plot([dspHRFprobe' ndspHRFprobe']); title('Constant onset, dispersion')
saveas(gcf,[opts.figdir,'HRFprobes.fig']);  
HRFprobe = [ onHRFprobe' dspHRFprobe' nonHRFprobe' ndspHRFprobe' ]; HRFprobe(isnan(HRFprobe)) = [];

%% use specified regressors to generate GS % residual BOLD
[cleanBOLD,res_ts,resBOLD] = genGS(sourceBOLD, mWBmask, M, HRFprobe, opts); %generate global signal regressor for removing global nuisance signals
if opts.denoise; cleanBOLD = denoiseData(cleanBOLD,WBmask,opts); end
if opts.smoothTS; cleanBOLD = smthData(cleanBOLD, mWBmask, opts); end

%% generate CVRindex map
opts.CVRidxdir = [opts.resultsdir,'CVRidx\']; mkdir(opts.CVRidxdir);
opts.fpass = [0.0001 0.08]; % https://doi.org/10.1148/radiol.2021203568
%opts.fpass = [0.2 0.3]; % https://doi.org/10.1177/0271678X20978582
[~, BP_ref, bpBOLD] = glmCVRidx(cleanBOLD, mWBmask, CERmask, opts);
save([opts.resultsdir,'processing_options.mat'], 'opts');

%end %for loop over subjects
