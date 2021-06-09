clear all
close all
global opts
    subj = 'sub_resp08'
    dir = ['D:\BACKUP\DATA\MRI_data\RO1\',subj,'\ses-7T\func\GE_TASK_RESP_5+10\'];
    cd(dir)
    seqpath = ['D:\BACKUP\DATA\MRI_data\RO1\',subj,'\ses-7T\derivatives\resp\'];
     
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
    maskBIN(maskBIN < 8) = 0; maskBIN(maskBIN >14) = 0;
    maskBIN(maskBIN > 0) = 1;
    
    
%directory & load options
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

    opts.TR = 0.85
    opts.seqpath = uigetdir(seqpath)
    opts.seqpath =  [opts.seqpath,'\'];
    [corrvec_CO2,corrvec_O2] = loadRAMRgen3(opts);
    % select hypercapnic block from CO2 trace
    figure; plot(corrvec_CO2/mean(corrvec_CO2),'r'); hold on; plot(corrvec_O2/mean(corrvec_O2),'b')
    title('isolate for timeseries correlation')
    [opts.traceidx,opts.traceidy] = ginput(2); opts.traceidx = round(opts.traceidx);
    close;
    CO2_corr = corrvec_CO2(1,opts.traceidx(1):opts.traceidx(2));
    O2_corr = corrvec_O2(1,opts.traceidx(1):opts.traceidx(2));

%% Align gas traces with BOLD data
cd(dir)
%generate timeseries for correlations
TS = meanTimeseries(BOLD,WBmask);
%GUI to correct alignment
trAlign(CO2_corr,O2_corr, TS,opts);
CO2trace = probe1;
O2trace = probe2;

%% Process Data
cd(dir)
opts.prefix = 'CO2BLOCK'
% setup function options
 
opts.denoise = 1; %perform wavelet denoising in the temporal dimension
opts.smoothTS = 1; %spatially smooth timeseries data
opts.smoothmap = 1; %option to smooth generated maps when needed for specific functions
opts.spatialdim = 2;% 2 for 2D smoothing 3 for 3D. Edge effects for 3D are handled with dilation algorithm
opts.FWHM = 2; %FWHM of smoothing filter

opts.motioncorr = 0.6; %motion parameters with higher correlation than threshold will not be included in the regression
opts.legOrder = [0]; %order of Legendre polynomials to include in regression (typically up to 4th order however can cause problems)
opts.plot = 0;


% pre-process data and run preliminary analysis

    opts.useGlobal = 0 %use global signal regressor during data-cleanup
    opts.priority = 1  %1 for CO2, 0 for O2
    useNuisance = 1; %use motion regressors  & legendre polynomials during data-cleanup
    
    cd(dir)
    if opts.priority
        probe = CO2trace;
        useNuisGas = 1; %use alternative gas as nuissance (most important when looking at O2 as primary stimulus
        nuiss_probe = O2trace;
    else
        probe = O2trace;
        useNuisGas = 1;
        nuiss_probe = CO2trace;
    end
    

    %update figure and results directories
    opts.resultsdir = [opts.savedir,opts.prefix,'_GS',int2str(opts.useGlobal),'_NS',int2str(useNuisance),'_WD',int2str(opts.denoise),'_CO2pr',int2str(opts.priority),'\'];
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
    mpfilename = ls('*-mp'); 
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
    opts.LVpercentile = 92
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
    autoCorr = abs(corr(meanTimeseries(rBOLD,maskBIN)',motionDer)); %remove
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
            np = []; np = [motionDer L nuiss_probe'];
        else
            np = []; np = [motionDer L];
            
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
    cleanBOLD = normTimeseries(cleanBOLD, mWBmask, [20 40]);
 
    %lag analysis 
    opts.trace_corr = 1;
    opts.corr_model = 1;
    opts.cvr_maps = 1;
    opts.glm_model = 1;
    opts.interp_factor = 5;
    opts.load_probe = 0;
    opts.highlag = 10;
    opts.upperlagthresh = 3;
    lagCVR(maskBIN.*mWBmask,mWBmask,cleanBOLD,nprobe,nuis_res, opts); %scrubbed data
    %lag analysis, GLM only using sourceBOLD and all nuissance regressors
    %lagCVR(GMmask.*mWBmask,mWBmask,cleanBOLD,nprobe,np, opts); %scrubbed data
    %save options
    save([opts.resultsdir,'processing_options.mat'], 'opts');
    


%% load lag map for sorting of clean data and production of carpet plots
%load lag map
[x y z dyn] = size(rBOLD);
lagdir = ['D:\BACKUP\DATA\MRI_data\RO1\sub_resp09\ses-7T\func\GE_TASK_RESP_5+10\CO2BLOCK_GS0_NS1_WD1_CO2pr1\corrLAG\']; cd(lagdir)
lagfile = 'robustLAG_r.nii.gz';
[lag_map,~,~] = loadImageData(lagdir, lagfile);
%load R2 map
r2file = 'r_map.nii.gz';
[r2_map,~,~] = loadImageData(lagdir, r2file);
r2_map(r2_map < 0.5) = 0; r2_map(r2_map > 0) = 1;
%load layer map
layerfile = 'varea.nii';
[varea,~,~] = loadImageData(dir, layerfile);
varea(varea > 1) = 0; varea(varea > 0) = 1;
mask = varea.*r2_map.*mWBmask;
masked_lag = lag_map.*mask; lag_mask = masked_lag; lag_mask(lag_mask~=0) = 1;
figure; subplot(1,2,1); imagesc(masked_lag(:,:,4));  subplot(1,2,2); imagesc(lag_mask(:,:,4));
masked_lag = masked_lag(:);
%masked_lag(mask == 0) = [];

[B,I] = sort(masked_lag);

%[voxels,coords] = grabTimeseries(rBOLD,lag_mask);
shag = reshape(rBOLD, [x*y*z dyn]);
shag = shag(I,:);
shag(B==0) = [];

%STD = std(svoxels,2);
svoxels = (shag-MEAN); 
carpetEdgeDetect(shag, opts.TR, 20,1,  .5,  0.40)