clear all
close all

%7t

dir = ['D:\BACKUP\DATA\MRI_data\CREATIV\BOLDmb\'];
seqpath  = ['D:\BACKUP\DATA\MRI_data\CREATIV\7tblock2\2021-05-06 09;58;52 (7Tramr)\'];


cd(dir)
%load data and masks
%filename = ls('*BOLD_masked_mcf.nii*')
filename = ls('*BOLD_applytopup.nii*')
[BOLD INFO BOLDheader] = loadImageData(dir, filename);

%file = ls('*seg1to*'); %GM segmentation
file = ls('*seg_0*'); %GM segmentation
[GMmask INFOmask HEADERmask] = loadImageData(dir, file);
GMmask(isinf(GMmask)) = 0;
file = ls('*seg_1*'); %WM segmentation
%file = ls('*seg2to*');%WM segmentation
[WMmask INFOmask HEADERmask] = loadImageData(dir, file);
WMmask(isinf(WMmask)) = 0;
file = ls('*seg_2*'); %CSF segmentation
%file = ls('*seg0to*');%CSF segmentation
[CSFmask INFOmask HEADERmask] = loadImageData(dir, file);
CSFmask(isinf(CSFmask)) = 0;

%for cases where you have to create your own mask
file = ls('*mean_brain_mask*') %WB segmentation
[WBmask INFOmask HEADERmask] = loadImageData(dir, file);
WBmask(isinf(WBmask)) = 0; WBmask = double(WBmask);

%directory & load options
global opts;

opts.TR = BOLDheader.dime.pixdim(5);
opts.voxelsize = BOLDheader.dime.pixdim(2:4);
opts.dyn = size(BOLD,4);
opts.xdata = [opts.TR:opts.TR:opts.TR*opts.dyn];

opts.savedir = dir;
opts.figdir = [dir,'figures\']; mkdir(opts.figdir);
opts.headers = struct();
opts.headers.ts = BOLDheader;
opts.headers.map = opts.headers.ts;
opts.headers.map.dime.dim(5) = 1;  
opts.headers.mask = HEADERmask;

%% load breathing traces

opts.seqpath = uigetdir(seqpath)

opts.remOUT = 1;
opts.remOUTmethod = 'median'; %'median' or 'mean' or 'quartiles' or 'grubbs' or 'gesd'
opts.remOUTbh = 0;
opts.evenstart = 'SequenceStart';
opts.eventend = 'SequenceEnd';
[CO2_corr,O2_corr] = loadRAMRgen4(opts); %loads RAMR gen4 gas traces
 
%% Align gas traces with BOLD data
cd(dir)
%generate timeseries for correlations
TS = meanTimeseries(BOLD,GMmask);

%GUI to correct alignment

[imf,residual,info] = vmd(TS)
trAlign(CO2_corr, O2_corr, mean(imf,2), opts);
CO2trace = probe1;
O2trace = probe2;

%test VMD




%% Process Data
cd(dir)
opts.prefix = 'CO2BLOCK' %prefix for save directory

%preprocessing options
opts.denoise = 1; %perform wavelet denoising in the temporal dimension
opts.smoothTS = 1; %spatially smooth timeseries data
opts.smoothmap = 1; %option to smooth generated maps when needed for specific functions
opts.useGlobal = 1 %use global signal regressor during data-cleanup
opts.priority = 1  %1 for CO2, 0 for O2
opts.motioncorr = 0.6; %motion parameters with higher correlation than threshold will not be included in the regression
opts.legOrder = [1 3]; %order of Legendre polynomials to include in regression (typically up to 4th order however can cause problems)

%smoothing options
opts.filtWidth = 5;
opts.spatialdim = 3;% 2 for 2D smoothing 3 for 3D. Edge effects for 3D are handled with dilation algorithm
opts.FWHM = 4; %FWHM of smoothing filter

opts.plot = 0;


% pre-process data and run preliminary analysis


    
    cd(dir)
    if opts.priority
        probe = CO2trace;
        useNuisGas = 0; %use alternative gas as nuissance (most important when looking at O2 as primary stimulus
        nuiss_probe = O2trace;
    else
        probe = O2trace;
        useNuisGas = 1;
        nuiss_probe = CO2trace;
    end
    

    %update figure and results directories
    opts.resultsdir = [opts.savedir,opts.prefix,'_GS',int2str(opts.useGlobal),'_WD',int2str(opts.denoise),'_CO2pr',int2str(opts.priority),'\'];
    mkdir(opts.resultsdir)
    opts.figdir = [opts.resultsdir,'figures\']; mkdir(opts.figdir);
    
    %reduce BOLD series to conserve memory for processing later on:
    %3rd argument can be start/end indices
    [idx, rBOLD] = chopTimeseries(BOLD,WBmask);

    %update options variables based on new epoch
    nuiss_probe = rescale(nuiss_probe(1,idx(1):idx(2)));
    nprobe = []; nprobe = probe(1,idx(1):idx(2));
    CO2_probe = CO2trace(1,idx(1):idx(2));
    O2_probe = O2trace(1,idx(1):idx(2));
    
    %account for signal onset
    opts.onset = 1:1:32;
    opts.disp = 1;
    opts.under = 1;
    [HRF,onHRFprobe] = convHRF(nprobe,opts); %Generate a double gamma HRF
    
    normprobe = onHRFprobe(1,:);
    onHRFprobe(1,:) = [];
    customMap = plasma(size(onHRFprobe,1));
    cshift = -round(5/opts.TR);
    conHRFprobe = [];
    for ii=1:size(onHRFprobe,1); conHRFprobe(ii,:) = circshift(onHRFprobe(ii,:),cshift); conHRFprobe(ii,end+cshift:end) = conHRFprobe(ii,end+cshift); end
    figure(30); subplot(1,2,1); hold on; for ii=1:size(customMap,1); plot(onHRFprobe(ii,:)', 'Color', customMap(ii,:)); end; plot(normprobe', 'k', 'LineWidth', 2); title('Shifted onset, limited dispersion')
    
    opts.onset = 1;
    opts.disp = 1:1:32
    opts.under = 1;
    [HRF,dspHRFprobe] = convHRF(nprobe, opts); %Generate a double gamma HRF
    
    %account for signal dispersion
    dspHRFprobe(1,:) = [];
    customMap = plasma(size(dspHRFprobe,1));
    figure(30); subplot(1,2,2); hold on; for ii=1:size(customMap,1); plot(dspHRFprobe(ii,:)', 'Color', customMap(ii,:)); end; plot(normprobe', 'k', 'LineWidth', 2); title('Dispersion')
    
    HRFprobe = [conHRFprobe' onHRFprobe' dspHRFprobe' ]
    [HRFprobe,~]=licols(HRFprobe); %remove linearly dependent components
    
    
    %load motion parameters
    mpfilename = ls('*mcf.par'); 
    nuisance = load(mpfilename); nuisance = nuisance(idx(1):idx(2),:);
    dtnuisance =  gradient(nuisance);
    %sqnuisance = nuisance.*nuisance;
    motionDer =[nuisance dtnuisance];

    %initialize Legendre regressors
    if opts.legOrder == 0; L = []; else
        L = []; for ii=opts.legOrder; L(ii,:) = rescale(LegendreN(ii,opts.xdata)); end; L = L';
    end
    %setup nuisance regressors to be used with HRFprobe to
    %determine global signal residual
    if useNuisGas
        np = []; np = [motionDer L nuiss_probe'];
    else
        np = []; np = [motionDer L];
    end
    %[np,~]=licols(np);
    
    % remove large vessel contributions that can affect CVR weigh
    % optimized BOLD regressor
    [mWBmask] = remLV(rBOLD,WBmask,opts);
    % remove signals explained by all known regressors to isolate
    % global signal fluctuations ('pseudo-resting state')
    [~,nuis_res,~] = genGS(rBOLD, mWBmask, np, HRFprobe, opts);
    %nuis_res = rescale(nuis_res,-1, 1);
    % Low frequency analysis on residual data
    %opts.fpass = [0.01 0.08]; %frequency range for ALFF/fALFF
   
    %fALFF(rBOLD, mWBmask, mWBmask, opts);
    % Generate CVR index map based on pseudo-RS or original data
    %opts.CVRidxdir = [opts.resultsdir,'CVRidx\']; mkdir(opts.CVRidxdir);
    % https://doi.org/10.1148/radiol.2021203568
    % https://journals.sagepub.com/doi/full/10.1177/0271678X16631755
    % https://doi.org/10.1177/0271678X20978582
    %opts.fpass = [0.001 0.1164];
    %[~, BP_WB, bpBOLD] = glmCVRidx(cleanBOLD, mWBmask,GMmask, opts);
    

    % Final cleanup of BOLD data using specified nuissance regressors
    %remove nuisance correlating with probe
    autoCorr = abs(corr(meanTimeseries(rBOLD,GMmask)',motionDer)); %remove
    autoCorr(autoCorr < opts.motioncorr) = 0; autoCorr(autoCorr > 0) = 1; %removes anything with more than weak correlation
    %remove highly correlated nuisance regressors to preserve signal response
    index = ([1:1:size(motionDer,2)]).*autoCorr; index(index == 0) = [];
    motionDer(:,index) = [];
    %motionDer = detrend(motionDer);   
    %setup nuisance regressors
    if opts.useGlobal
        if useNuisGas
            np = []; np = [motionDer L nuis_res' nuiss_probe'];
        else
            np = []; np = [motionDer L nuis_res'];
        end
    else
        if useNuisGas
            np = []; np = [motionDer L nuis_res' nuiss_probe'];
        else
            np = []; np = [motionDer L nuis_res'];
            
        end
    end
    %[np,~]=licols(np);
  
    %remove signals explained by nuissance regressors
    [cleanBOLD] = scrubData(rBOLD, mWBmask, np,nprobe, opts);
    %denoise
    if opts.denoise; [cleanBOLD] = denoiseData(cleanBOLD,WBmask,opts); end
    %spatial smoothing
    if opts.smoothTS; cleanBOLD = smthData( cleanBOLD, mWBmask, opts); end
    %normalize
    cleanBOLD = normTimeseries(cleanBOLD, mWBmask);
    %save cleaned timesieres
    %lag analysis 
    opts.trace_corr = 1;
    opts.corr_model = 1;
    opts.cvr_maps = 1;
    opts.glm_model = 1;
    opts.interp_factor = 5;
    opts.lowerlagthresh = -1
    opts.opperlagthresh = 3
    opts.lowlag = -3
    opts.highlag = 45
    lagCVR(GMmask.*mWBmask,mWBmask,cleanBOLD,nprobe,nuiss_probe, opts); %scrubbed data
    %lag analysis, GLM only using sourceBOLD and all nuissance regressors
    %lagCVR(GMmask.*mWBmask,mWBmask,cleanBOLD,nprobe,np, opts); %scrubbed data
    %save options
    save([opts.resultsdir,'processing_options.mat'], 'opts');
    
