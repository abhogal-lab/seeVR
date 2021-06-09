clear all
close all
RO1 = 1;
if RO1
    
    
    subj = 'sub_resp09'
    dir = ['D:\BACKUP\DATA\MRI_data\RO1\',subj,'\ses-7T\func\GE_TASK_RESP_5+10\'];
    cd(dir)
    seqpath = ['D:\BACKUP\DATA\MRI_data\RO1\',subj,'\ses-7T\derivatives\resp\'];
    GEN4=0;
    
    
    %load data and masks
    filename = ls('*-mc-w.nii*')
    [BOLD INFO BOLDheader] = loadImageData(dir, filename);
    BOLD(BOLD == 0) = NaN; BOLD(isinf(BOLD)) = NaN; BOLD = double(BOLD);
    
    file = ls('*-mask.nii*'); %segmentation
    [WBmask INFOmask HEADERmask] = loadImageData(dir, file);
    
    file = ls('*lay.nii*'); %segmentation
    [laymask INFOmask HEADERmask] = loadImageData(dir, file);
    mask = laymask;
    mask(isinf(mask)) = 0;
    maskBIN = mask;
    maskBIN(maskBIN < 14) = 0;
    maskBIN(maskBIN > 0) = 1;
    
    GMmask = maskBIN;
    %generate timeseries for correlations
    TS = make_mean_ts(single(BOLD),maskBIN);
    
else
    
    %activate other types of data
    %CHAMPAGNE
    
    dir = ['C:\Users\abhogal\Documents\DATA\MRI_data\Champagne\Subj2\bh.exp.x4.bold.mb.exp\'];
    seqpath  = ['C:\Users\abhogal\Documents\DATA\MRI_data\Champagne\Subj2\ramr\'];
    
    
    
    %APRICOT
    %     vnumber = ['APP003']
    %     dir = ['C:\Users\abhogal\Documents\DATA\MRI_data\APRICOT\',vnumber,'\BOLD\'];
    %     seqpath  = ['C:\Users\abhogal\Documents\DATA\MRI_data\APRICOT\',vnumber,'\traces\'];
    %     opts.GEN4 = 1;
    
    
    
    %7T data
    %         vn = 'VR001'
    %         vn(vn == ' ') = [];
    %         type = 'BLOCK';
    %         dir = ['C:\Users\abhogal\Documents\DATA\MRI_data\delayCOMP\',vn,'\',type,'\'];
    %         cd(dir)
    %         seqpath = ['C:\Users\abhogal\Documents\DATA\MRI_data\delayCOMP\',vn,'\Sequences\'];
    %         GEN4 = 0;
    %
    %
    %     OTHER (i.e. SINUS)
    %     dir = 'D:\BACKUP\DATA\MRI_data\SINUS\V3392\ANALYSISCO2\';
    %     seqpath = 'D:\BACKUP\DATA\MRI_data\SINUS\V3392\Sequences\';
    %     GEN4 = 0;
    %
    
    makeWBmask = 0;
    
    cd(dir)
    %load data and masks
    filename = ls('*BOLD_masked_mcf.nii*')
    %filename = ls('*BOLD_applytopup.nii*')
    [BOLD INFO BOLDheader] = loadImageData(dir, filename);
    
    file = ls('*seg_0*'); %GM segmentation
    [GMmask INFOmask HEADERmask] = loadImageData(dir, file);
    GMmask(isinf(GMmask)) = 0;
    file = ls('*seg_1*');%WM segmentation
    [WMmask INFOmask HEADERmask] = loadImageData(dir, file);
    WMmask(isinf(WMmask)) = 0;
    
    file = ls('*seg_2*');%CSF segmentation
    [CSFmask INFOmask HEADERmask] = loadImageData(dir, file);
    CSFmask(isinf(CSFmask)) = 0;
    
    %for cases where you have to create your own mask
    if makeWBmask
        tmp = mean(BOLD,4);
        [WBmask] = make_mask(mean(BOLD,4)); clear tmp
        saveImageData(WBmask, HEADERmask, savedir,  'BOLD_mean_brain_mask.nii.gz', INFOmask.dim, INFOmask.datatype);
    else
        file = ls('*mean_brain_mask*') %WB segmentation
        [WBmask INFOmask HEADERmask] = loadImageData(dir, file);
        WBmask(isinf(WBmask)) = 0; WBmask = double(WBmask);
    end
    
end

%directory & load options
if RO1; opts.GEN4 = 0; else opts.GEN4 = 1; end
opts.TR = BOLDheader.dime.pixdim(5);
opts.voxelsize = BOLDheader.dime.pixdim(2:4);
opts.dyn = size(BOLD,4); xdata = [opts.TR:opts.TR:opts.TR*opts.dyn];
opts.savedir = [dir];
opts.trfigdir = [dir,'tr_align\']; mkdir(opts.trfigdir);
opts.headers = struct();
opts.headers.ts = BOLDheader;
opts.headers.map = HEADERmask;
%% load breathing traces

if opts.GEN4
    opts.seqpath = uigetdir(seqpath)
    opts.extra = 30; %defines number of extra samples taken before and after the start and end of the events
    opts.remOUT = 1;
    opts.remOUTbh = 1;
    opts.evenstart = 'SequenceStart';
    opts.eventend = 'SequenceEnd';
    [corrvec_CO2,corrvec_O2] = loadRAMRgen4(opts); %loads RAMR gen4 gas traces
    CO2_corr = corrvec_CO2;
    O2_corr = corrvec_O2;
    x_indyi(1) = 1; x_indyi(2)= length(CO2_corr);
else
    opts.TR = 0.8
    opts.seqpath = uigetdir(seqpath)
    opts.seqpath =  [opts.seqpath,'\'];
    [corrvec_CO2,corrvec_O2] = loadRAMRgen3(opts);
    % select hypercapnic block from CO2 trace
    figure; plot(corrvec_CO2/mean(corrvec_CO2),'r'); hold on; plot(corrvec_O2/mean(corrvec_O2),'b')
    title('isolate for timeseries correlation')
    [x_indyi,y_ind] = ginput(2); x_indyi = round(x_indyi);
    close;
    CO2_corr = corrvec_CO2(1,x_indyi(1):x_indyi(2));
    O2_corr = corrvec_O2(1,x_indyi(1):x_indyi(2));
end

%% Align gas traces with BOLD data
cd(dir)
%generate timeseries for correlations
TS = make_mean_ts(BOLD,GMmask);
%GUI to correct alignment
trAlign(corrvec_CO2, corrvec_O2, TS, x_indyi, opts.trfigdir);

%% Process Data
cd(dir)
% setup function options
opts.prewhite = 0; %zero mean and unit variance of data. normally leads to bad results
opts.interp_factor = 4; %factor by which to temporally interpolate data. Better for picking up lags between TR
opts.denoise = 1; %perform wavelet denoising in the temporal dimension
opts.smoothTS = 1; %spatially smooth timeseries data
opts.spatialdim = 2;% 2 for 2D smoothing 3 for 3D. Edge effects for 3D are handled with dilation algorithm
opts.FWHM = 1.5; %FWHM of smoothing filter
opts.save_rts = 0; %save correlation timeseries - can be used to visualize lags
opts.corr_model = 1; %perform correlation based analysis
opts.glm_model = 1; %perform GLM based lag regression. This works well for high temporal resolution data
opts.cvr_maps = 1; %generate CVR maps based on regression of BOLD signal with gas probe
opts.effective_ET = 1; %uses effective CO2 to calculate CVR
opts.uni = 1;% if this is set to 1 negative correlations will be ignored
opts.loadprobe = 0; %saves time for creating regressor if one with the correct length already exists
opts.TraceCorr = 1; %perform additional correlation with gas trace (on top of optimized regressor)
opts.LVthresh = 0.18; %threshold for removing bright signals via updated WB mask; this parameters is dataset-specific
opts.motioncorr = 0.4; %motion parameters with higher correlation than threshold will not be included in the regression
opts.legOrder = [1 3]; %order of Legendre polynomials to include in regression (typically up to 4th order however can cause problems)
opts.plot = 0;
opts.saveCleaned = 0; %set to 1 to save BOLD data with motion/legendre regressed out
%by setting opts.refine_regressor = 0 with opts.TraceCorr = 1 analysis will
%be done on end-tidal traces only

%%%refinement optimized regressor
opts.refine_regressor = 1; %refine BOLD regressor. If '0' a straight correlation will be done with probe
opts.pca_analysis= 1; %PCA analysis to optimize BOLD regressor
opts.corrthresh = 0.4; %threshold by which to accept correlated voxels during the refinement of the BOLD regressor
opts.norm_regr = 0; %normalize the regressor for correlation analysis between 0-1
%thresholds for refinening optimized regressors (+/- 3 has given good
%results!)
opts.lowerlagthresh = -2;
opts.upperlagthresh = 2;
%thresholds for lag map creation
opts.lowlag = -8; %setup lower lag limit; negative for misalignment and noisy correlation
opts.highlag = 75; %setups upper lag limit; allow for long lags associated with pathology

opts.lowlag = round((opts.lowlag/opts.TR)*opts.interp_factor); %setup lower lag limit; negative for misalignment and noisy correlation
opts.highlag = round((opts.highlag/opts.TR)*opts.interp_factor); %setups upper lag limit; allow for long lags associated with pathology
opts.lowerlagthresh = round(opts.lowerlagthresh/opts.TR*opts.interp_factor);
opts.upperlagthresh = round(opts.upperlagthresh/opts.TR*opts.interp_factor);

%% pre-process data and run preliminary analysis

for useGlobal=[1 0] %use global signal regressor during data-cleanup
    opts.priority= 1  %1 for CO2, 0 for O2
    useNuisance = 1; %use motion regressors  & legendre polynomials during data-cleanup
    
    cd(dir)
    if priority
        probe = CO2trace;
        useNuisGas = 1; %use alternative gas as nuissance (most important when looking at O2 as primary stimulus
        nuiss_probe = O2trace;
        prefix = '4CO2BLOCK'
    else
        probe = O2trace;
        useNuisGas = 1;
        nuiss_probe = CO2trace;
        prefix = 'O2BLOCK'
        
    end
    
    %update figure and results directories
    opts.resultsdir = [opts.savedir,prefix,'_GS',int2str(useGlobal),'_NS',int2str(useNuisance),'_WD',int2str(denoise),'_CO2pr',int2str(priority),'\'];
    mkdir(opts.resultsdir)
    opts.figdir = [opts.resultsdir,'figures\']; mkdir(opts.figdir);
    
    %reduce BOLD series to conserve memory for processing later on:
    %3rd argument can be start/end indices
    [idx, rBOLD] = chop_ts(BOLD,WBmask,[100 580]);
    %update options variables based on new epoch
    opts.dyn = []; [opts.xdim,opts.ydim,opts.zdim,opts.dyn] = size(rBOLD)
    xdata = []; xdata = [opts.TR:opts.TR:opts.TR*opts.dyn];
    nuiss_probe = rescale(nuiss_probe(1,idx(1):idx(2)));
    nprobe = []; nprobe = probe(1,idx(1):idx(2));
    CO2_probe = CO2trace(1,idx(1):idx(2));
    O2_probe = O2trace(1,idx(1):idx(2));
    
    %account for signal onset
    opts.rratio = 1000; %seems to affect the vertical spread of HRF (use large value to limit)
    opts.onset = 1:1*30;
    opts.disp = 0.4;
    opts.under = 1;
    [HRF,onHRFprobe] = convHRF(nprobe,opts); %Generate a double gamma HRF
    
    normprobe = onHRFprobe(1,:);
    onHRFprobe(1,:) = [];
    cshift = -round(5/opts.TR);
    for ii=1:size(onHRFprobe,1); onHRFprobe(ii,:) = circshift(onHRFprobe(ii,:),cshift); onHRFprobe(ii,end+cshift:end) = onHRFprobe(ii,end+cshift); end
    figure(30); subplot(1,2,1); plot(onHRFprobe'); hold on; plot(normprobe', 'k', 'LineWidth', 2); title('Shifted onset, limited dispersion')
    
    opts.onset = 1;
    opts.disp = 1:1*20
    opts.under = 1;
    [HRF,dspHRFprobe] = convHRF(nprobe, opts); %Generate a double gamma HRF
    
    %account for signal dispersion
    dspHRFprobe(1,:) = [];
    figure(30); subplot(1,2,2); plot(dspHRFprobe'); hold on; plot(normprobe', 'k', 'LineWidth', 2); title('Dispersion')
    
    HRFprobe = [ onHRFprobe' dspHRFprobe' ]
    [HRFprobe,~]=licols(HRFprobe); %remove linearly dependent components
    
    
    %load motion parameters
    if RO1; mpfilename = ls('*-mp'); else; mpfilename = ls('*mcf.par'); end
    nuisance = load(mpfilename); nuisance = nuisance(idx(1):idx(2),:);
    dtnuisance =  gradient(nuisance);
    sqnuisance = nuisance.*nuisance;
    motionDer =[nuisance dtnuisance sqnuisance];
    
    %initialize Legendre regressors
    if opts.legOrder == 0; L = []; else
        L = []; for ii=opts.legOrder; L(ii,:) = rescale(LegendreN(ii,xdata)); end
    end
    %setup nuisance regressors to be used with HRFprobe to
    %determine global signal residual
    if useNuisGas
        np = []; np = [motionDer L' nuiss_probe'];
    else
        np = []; np = [motionDer L'];
    end
    [np,~]=licols(np);
    
    % remove large vessel contributions that can affect CVR weigh
    % optimized BOLD regressor
    [mWBmask] = remLV(rBOLD,WBmask,opts);
    % remove signals explained by all known regressors to isolate
    % global signal fluctuations ('pseudo-resting state')
    [~,nuis_res,resBOLD] = genGS(rBOLD, mWBmask, np,HRFprobe, opts);
    % Low frequency analysis on residual data
    opts.fpass = [0.01 0.08]; %frequency range for ALFF/fALFF
    opts.ALFFdir = [opts.resultsdir,'ALFF\']; mkdir(opts.ALFFdir);
    fALFF(rBOLD, mWBmask, mWBmask, opts);
    % Generate CVR index map based on pseudo-RS or original data
    opts.CVRidxdir = [opts.resultsdir,'CVRidx\']; mkdir(opts.CVRidxdir);
    % https://doi.org/10.1148/radiol.2021203568
    % https://journals.sagepub.com/doi/full/10.1177/0271678X16631755
    % https://doi.org/10.1177/0271678X20978582
    opts.fpass = [0.001 0.1164];
    [~, BP_WB, bpBOLD] = glmCVRidx(resBOLD, mWBmask,mWBmask, opts);
    
    clear resBOLD
    
    % Final cleanup of BOLD data using specified nuissance regressors
    %remove nuisance correlating with probe
    autoCorr = abs(corr(make_mean_ts(rBOLD,mWBmask)',motionDer)); %remove
    autoCorr(autoCorr < opts.motioncorr) = 0; autoCorr(autoCorr > 0) = 1; %removes anything with more than weak correlation
    %remove highly correlated nuisance regressors to preserve signal response
    index = ([1:1:size(motionDer,2)]).*autoCorr; index(index == 0) = [];
    motionDer(:,index) = [];
    
    %setup nuisance regressors
    if useGlobal
        if useNuisGas
            np = []; np = [motionDer L' nuis_res' nuiss_probe'];
        else
            np = []; np = [motionDer L' nuis_res'];
        end
    else
        if useNuisGas
            np = []; np = [L' nuis_res' nuiss_probe'];
        else
            np = []; np = [L' nuis_res'];
            
        end
    end
    [np,~]=licols(np);
    
    %remove signals explained by nuissance regressors
    [cleanBOLD] = scrubBOLD(rBOLD, mWBmask, nprobe', np, opts.figdir);
    
    if opts.denoise; [cleanBOLD] = denoiseData(cleanBOLD,WBmask,opts); end
    
    
    if opts.smoothTS
        cleanBOLD( cleanBOLD == 0) = NaN;
        cleanBOLD = smthData( cleanBOLD, mWBmask, opts);
    end
    %normalize processed timeseries to baseline (if no index
    %supplied in third argument, user may choose begin/end baseline
    cleanBOLD = norm_ts(cleanBOLD, mWBmask, [5 30]);
    %save cleaned timesieres
    if opts.saveCleaned; saveImageData(cleanBOLD, BOLDheader, opts.resultsdir, 'cleanBOLD.nii.gz', 64); end
    %setup lag analysis directories
    lagCVR(GMmask.*mWBmask,mWBmask,cleanBOLD,nprobe,nuiss_probe', opts); %scrubbed data
    %lagCVR(GMmask.*mWBmask,mWBmask,cleanBOLD,nprobe,np, opts); %for GLM only add nuissance regressors
    %clear voxel_ts error_data np_coef nuisance dtnuisance cleanBOLD nnuiss nnnuiss nnuis_TS nnnuis_TS probe nuis_probe
    save([opts.resultsdir,'processing_options.mat'], 'opts');
    
end

%% Generate CVR, lag, corrected CVR maps (correlation & GLM based)

loadprobe = 1; %use BOLD regressor generated above instead of endtidal trace
smoothTempHRF =1;
HRF_mapsRO1(TR,headers,GMmask.*mWBmask,mWBmask,BOLD_ts,probe(1,idx(1):idx(2))',smoothTempHRF,prewhite,interp_factor,dim,savedir1,nuiss_probe',loadprobe)
load gong
sound(y,Fs)


%% Generate dynamic movie
funcopts.scale = [-5 5];
funcopts.row = 5;
funcopts.col = 6;
funcopts.step = 1;
funcopts.start = 14;
funcopts.start_ind = 24;
funcopts.step = 1;
funcopts.type = 0; %0 for no output, 1 for .gif, 2 for .mp4

generateMovie(funcopts,WBmask,BOLD,imrotate(cleanBOLD,90))
