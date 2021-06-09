clear all
close all
RO1 = 0;
if RO1
    
   
    subj = 'sub_resp08'
    dir = ['D:\BACKUP\DATA\MRI_data\RO1\',subj,'\ses-7T\func\GE_TASK_RESP_5+10\'];
    cd(dir)
    seqpath = ['D:\BACKUP\DATA\MRI_data\RO1\',subj,'\ses-7T\derivatives\resp\'];
    GEN4=0;
    
    
    %load data and masks
    filename = ls('*-mc-w.nii*')
    [BOLD INFO BOLDheader] = loadImageData(dir, filename);
    BOLD(BOLD == 0) = NaN; BOLD(isinf(BOLD)) = NaN; BOLD = double(BOLD);
    dyn = size(BOLD,4); xdim = size(BOLD,1); ydim = size(BOLD,2); zdim = size(BOLD,3);
    TR = 0.85;
    xdata = [TR:TR:TR*dyn];
    dim = [1 1 1];
    
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
%     
%     dir = ['C:\Users\abhogal\Documents\DATA\MRI_data\Champagne\mri\bh.exp.x4.bold.mb.insp\'];
%     seqpath  = ['C:\Users\abhogal\Documents\DATA\MRI_data\Champagne\ramr\'];
%     GEN4 = 1;
    
    
    %APRICOT
    vnumber = ['APP007']
    dir = ['C:\Users\abhogal\Documents\DATA\MRI_data\APRICOT\',vnumber,'\ANALYSIS\'];
    seqpath  = ['C:\Users\abhogal\Documents\DATA\MRI_data\APRICOT\',vnumber,'\traces\'];
    GEN4 = 1;
    TR = 1.05;
    dim = [2.5 2.5 2.5];
    
    %7T data
    %         vn = 'VR001'
    %         vn(vn == ' ') = [];
    %         type = 'BLOCK';
    %         dir = ['C:\Users\abhogal\Documents\DATA\MRI_data\delayCOMP\',vn,'\',type,'\'];
    %         cd(dir)
    %         seqpath = ['C:\Users\abhogal\Documents\DATA\MRI_data\delayCOMP\',vn,'\Sequences\'];
    %         GEN4 = 0;
    %         TR = 3;
    %         dim = [1.5 1.5 1.6];
    %
    %     OTHER (i.e. SINUS)
    %     dir = 'D:\BACKUP\DATA\MRI_data\SINUS\V3392\ANALYSISCO2\';
    %     seqpath = 'D:\BACKUP\DATA\MRI_data\SINUS\V3392\Sequences\';
    %     GEN4 = 0;
    %     TR = 0.7;
    %     dim = [1.042  1.042  1.3];
    
    makeWBmask = 0;
        
    cd(dir)
    %load data and masks
    filename = ls('*masked_mcf.nii*')
    %filename = ls('BOLD_test_mut_refmean_st4_edge.nii') %for V2407
    [BOLD INFO BOLDheader] = loadImageData(dir, filename); BOLD = double(BOLD);
    opts.TR = BOLDheader.dime.pixdim(5);
    opts.dyn = size(BOLD,4); xdata = [opts.TR:opts.TR:opts.TR*opts.dyn];
    
    file = ls('*seg_0*'); %GM segmentation
    [GMmask INFOmask HEADERmask] = loadImageData(dir, file);
    GMmask(isinf(GMmask)) = 0; GMmask = double(GMmask);
    
    file = ls('*seg_1*');%WM segmentation
    [WMmask INFOmask HEADERmask] = loadImageData(dir, file);
    WMmask(isinf(WMmask)) = 0; WMmask = double(WMmask);
    
    file = ls('*seg_2*');%CSF segmentation
    [CSFmask INFOmask HEADERmask] = loadImageData(dir, file);
    CSFmask(isinf(CSFmask)) = 0; CSFmask = double(CSFmask);
    
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

opts.savedir = [dir];
opts.trfigdir = [dir,'tr_align\']; mkdir(opts.trfigdir);
opts.headers = struct();
opts.headers.ts = BOLDheader;
opts.headers.map = HEADERmask;
%% load breathing traces

if GEN4
    opts.seqpath = uigetdir(seqpath)
    cd(opts.seqpath)
    %%step 1: import breathing files
    
    
    %Import EndTidal
    file = ls('EndTidal.*')
    filename = [opts.seqpath,'/',file];
    [MRTimes,DesiredPO2mmHg,DesiredPCO2mmHg,AchievablePO2mmHg,AchievablePCO2mmHg,PO2mmHg,PCO2mmHg,RestingPO2mmHg,RestingPCO2mmHg,PBarommHg,Inspiretimeseconds,Expiretimeseconds,Breathidx,TidalvolumemL,RespirationrateBPM,StartInspiresec,O2AdjustmentmmHg,CO2AdjustmentmmHg,G1TargetvolmL,G1FCO2,G1FO2,G2FCO2,G2FO2] = import_EndTidal(filename);
    
    %import RGM
    file = ls('RGM.*')
    filename = [opts.seqpath,'/',file];
    [MRTimes1,PO2mmHg1,PCO2mmHg1,PBarommHg1,PMouthmmH2O,FlowMouthmLmin,FlowS1mLmin,FlowS2mLmin,BreathPhase] = import_RGM(filename);
    
    %import events
    file = ls('Events.*')
    filename = [opts.seqpath,'/',file];
    [MRTimes3, CtrlRoomTimes, Event] = import_Events(filename);
    
    %import physiological parameters
    file = ls('PhysioParams.*')
    filename = [opts.seqpath,'/',file];
    [MRTimes2, ID, FRCmL, VdmL, TissuestoreO2mL, TissuestoreCO2mL, VO2mLmin, VCO2mLmin, QmLmin, hBconcentrationgdLBlood, Restingmetabolicscalefactor, ResponseReason] = import_Physio(filename)
    
    % resample and realign the breathing trace and the Endtidal trace to have the same start and end and same sampling rate
    opts.extra = 30; %defines number of extra samples taken before and after the start and end of the events
    opts.remOUT = 1;
    opts.remOUTbh = 0;
    opts.evenstart = 'SequenceStart';
    opts.eventend = 'SequenceEnd';
    opts.remOUTmethod = 'quartiles'; %'median' or 'mean' or 'quartiles' or 'grubbs' or 'gesd'
    [nxi,corrvec_CO2,corrvec_O2,nxi1,rawCO2,rawO2] = resampleEndtidalBreathing(MRTimes,PCO2mmHg,PO2mmHg,MRTimes1,PCO2mmHg1,PO2mmHg1,MRTimes3,Event,opts)  
    
else
    cd(opts.seqpath)
    
    oversample = 1;
    [BBBdata corrvec_CO2 corrvec_O2] = resampleGEN3(opts.seqpath,TR,dyn,oversample)
    %import RAW
    rawFile = ls('*raw*')
    raw = importRAW(rawFile);
    
    %plot time versus BBB, raw is: 1)time, 2)Pmouth, 3)Pdelay, 4)PCO2, 5)PO2
    figure;
    sz = 140;
    plot(table2array(raw(:,1)),table2array(raw(:,4)),'k', 'LineWidth', 0.5); hold on;
    scatter(BBBdata.data(:,1),BBBdata.data(:,2),sz, 'm.');
    figure; %O2
    sz = 140;
    plot(table2array(raw(:,1)),table2array(raw(:,5)),'k', 'LineWidth', 0.5); hold on
    scatter(BBBdata.data(:,1),BBBdata.data(:,3),sz, 'c.');
    
end

%%
if GEN4
    CO2_corr = corrvec_CO2;
    O2_corr = corrvec_O2;
    x_indyi(1) = 1; x_indyi(2)= length(CO2_corr);
else
    % select hypercapnic block from CO2 trace
    figure; plot(corrvec_CO2/mean(corrvec_CO2),'r'); hold on; plot(corrvec_O2/mean(corrvec_O2),'b')
    title('isolate for timeseries correlation')
    [x_indyi,y_ind] = ginput(2); x_indyi = round(x_indyi);
    close;
    CO2_corr = corrvec_CO2(1,x_indyi(1):x_indyi(2));
    O2_corr = corrvec_O2(1,x_indyi(1):x_indyi(2));
end
cd(dir)
% check alignment
%generate timeseries for correlations
TS = make_mean_ts(BOLD,GMmask);

%GUI to correct alignment
trAlign(corrvec_CO2, corrvec_O2, TS, x_indyi, opts.trfigdir);

%% initial data 'scrubbing'
cd(dir)
% setup function options
opts.prewhite = 0; %zero mean and unit variance of data. normally leads to bad results
opts.interp_factor = 4; %factor by which to temporally interpolate data. Better for picking up lags between TR
opts.spatialdim = 3;% 2 for 2D smoothing 3 for 3D. Edge effects for 3D are handled with dilation algorithm
opts.FWHM = 3.5; %FWHM of smoothing filter
opts.save_rts = 0; %save correlation timeseries - can be used to visualize lags
opts.pca_analysis=1; %PCA analysis to optimize BOLD regressor
opts.norm_regr = 0; %normalize the regressor for correlation analysis between 0-1
opts.glm_model = 1; %perform GLM based lag regression. This works well for high temporal resolution data
opts.cvr_maps = 1; %generate CVR maps based on regression of BOLD signal with gas probe
opts.effective_ET = 1; %uses effective CO2 scaled to represent BOLD regressor
opts.uni = 0;% if this is set to 1 negative correlations will be ignored
opts.loadprobe = 0; %saves time for creating regressor if one already exists
opts.TraceCorr = 1; %perform additional correlation with gas trace (on top of optimized regressor)
opts.corrthresh = 0.5; %threshold by which to accept correlated voxels during the refinement of the BOLD regressor
opts.LVthresh = 0.12; %threshold for removing bright signals via updated WB mask; this parameters is dataset-specific 
opts.motioncorr = 0.3; %motion parameters with higher correlation than threshold will not be included in the regression
opts.refine_regressor = 1; %refine BOLD regressor. If '0' a straight correlation will be done with probe
%% pre-process data and run preliminary analysis
for priority=[1] %priority = aa %1 for CO2, 0 for O2
    for denoise=[1] %perform wavelet denoising in the temporal dimension
        for useGlobal=[1] %use global signal regressor during data-cleanup
            
            smthdata = 1; %turn on spatial smoothing
            regrNuiss = 0; %regress out non-prioroty gas trace - mainly useful for strong CO2 effects on O2
            save_cleaned = 0; %set to 1 to save BOLD data with motion/legendre regressed out
            smoothMP = 0; %smooth motion parameters before regression
            useNuisGas = 0; %use alternative gas as nuissance (most important when looking at O2 as primary stimulus)
            useNuisance = 1; %use motion regressors  & legendre polynomials during data-cleanup
            
            cd(dir)
            if priority
                probe = CO2trace;
                nuiss_probe = O2trace;
                prefix = '4CO2BLOCK'
%                                  x_ind(1) = 10;%hard code or use the reduce code below
%                                  x_ind(2) = size(BOLD,4);
            else
                probe = O2trace;
                nuiss_probe = CO2trace;
                prefix = 'O2BLOCK'
                %x_ind(1) = 690;
                %x_ind(2) = size(BOLD,4);
            end
            opts.resultsdir = [opts.savedir,prefix,'_GS',int2str(useGlobal),'_NS',int2str(useNuisance),'_WD',int2str(denoise),'_CO2pr',int2str(priority),'\'];
            mkdir(opts.resultsdir)
            opts.figdir = [opts.resultsdir,'figures\']; mkdir(opts.figdir);
            
            %reduce BOLD series to conserve memory for processing later on
            figure;  plot(make_mean_ts( BOLD, WBmask));
            if priority; title('Select Epoch'); else;  title('Select O2 Epoch'); end
            [x_ind,y_ind] = ginput(2); x_ind = round(x_ind); %select points
            close;
            
            rBOLD = BOLD(:,:,:,x_ind(1):x_ind(2));
            opts.dyn = []; [xdim ydim zdim dyn] = size(rBOLD)
            
            xdata = []; xdata = [opts.TR:opts.TR:opts.TR*opts.dyn];
            
            nuiss_probe = rescale(nuiss_probe(1,x_ind(1):x_ind(2)));
            nprobe = []; nprobe = probe(1,x_ind(1):x_ind(2));
            CO2_probe = CO2trace(1,x_ind(1):x_ind(2)); %extract CO2 info
            O2_probe = O2trace(1,x_ind(1):x_ind(2)); %extract O2 info
            
            %account for signal onset             
            onset_vec = 1:1*30
            disp_vec = 0.4
            undersht_vec = 1
            [HRF,onHRFprobe] = convHRF(nprobe, onset_vec, disp_vec, undersht_vec); %Generate a double gamma HRF

            normprobe = onHRFprobe(1,:);
            onHRFprobe(1,:) = [];
            cshift = -round(5/opts.TR);
            for ii=1:size(onHRFprobe,1); onHRFprobe(ii,:) = circshift(onHRFprobe(ii,:),cshift); onHRFprobe(ii,end+cshift:end) = onHRFprobe(ii,end+cshift); end
            figure(30); subplot(1,2,1); plot(onHRFprobe'); hold on; plot(normprobe', 'k', 'LineWidth', 2); title('Shifted onset, limited dispersion')
            
            onset_vec = 1;
            disp_vec = 1:2*24
            undersht_vec = 1;
            [HRF,dspHRFprobe] = convHRF(nprobe, onset_vec, disp_vec, undersht_vec); %Generate a double gamma HRF

            %account for signal dispersion
            dspHRFprobe(1,:) = [];
            figure(30); subplot(1,2,2); plot(dspHRFprobe'); hold on; plot(normprobe', 'k', 'LineWidth', 2); title('Dispersion')
            
            HRFprobe = [ onHRFprobe' dspHRFprobe' ]
            [HRFprobe,idx]=licols(HRFprobe); %remove linearly dependent components
            
            
            %load motion parameters
            if RO1; mpfilename = ls('*-mp'); else; mpfilename = ls('*mcf.par'); end
            nuisance = load(mpfilename); nuisance = nuisance(x_ind(1):x_ind(2),:);
            dtnuisance = gradient(nuisance);
            sqnuisance = nuisance.*nuisance;
            %initialize Legendre regressors
            L = []; for ii=1:3; L(ii,:) = rescale(LegendreN(ii,xdata)); end
            
            if useNuisGas
                np = []; np = [nuisance dtnuisance sqnuisance L' nuiss_probe'];            
            else
                np = []; np = [nuisance dtnuisance sqnuisance L'];
            end
            [np,idx]=licols(np); %remove linearly dependent components
            
            % remove large vessel contributions that can affect CVR weigh
            % optimized BOLD regressor
            [mWBmask] = remLV(rBOLD,WBmask,opts);
            % remove signals explained by all known regressors to isolate
            % global signal fluctuations ('pseudo-resting state')
            [~,resBOLD] = genGS(rBOLD, mWBmask, np,[], opts);
            nuis_res = rescale(make_mean_ts(resBOLD, mWBmask));
            % Low frequency analysis on residual data
            opts.fpass = [0.01 0.08]; %frequency range for ALFF/fALFF
            opts.ALFFdir = [opts.resultsdir,'ALFF\']; mkdir(opts.ALFFdir);
            fALFF(rBOLD, mWBmask, mWBmask, opts);
            % Generate CVR index map based on pseudo-RS or original data
            opts.CVRidxdir = [opts.resultsdir,'CVRidx\']; mkdir(opts.CVRidxdir);
            %opts.fpass = [0.001 0.1164]; % https://doi.org/10.1148/radiol.2021203568
            opts.fpass = [0.001 0.1164]; % https://doi.org/10.1177/0271678X20978582
            [~, BP_WB, bpBOLD] = glmCVRidx(resBOLD, mWBmask,mWBmask, opts);
            
            clear resBOLD
            
            % Final cleanup of BOLD data using specified nuissance regressors
            %remove nuisance correlating with probe
            autoCorr = abs(corr(make_mean_ts(rBOLD,mWBmask)',nuisance)); %remove
            autoCorr(autoCorr < opts.motioncorr) = 0; autoCorr(autoCorr > 0) = 1; %removes anything with more than weak correlation
            %remove highly correlated nuisance regressors to preserve signal response
            index = ([1:1:size(nuisance,2)]).*autoCorr; index(index == 0) = [];
            nuisance(:,index) = []; dtnuisance(:,index) = []; sqnuisance(:,index) = [];
            
            if smoothMP %smooth motion params
                for ii=1:size(nuisance,2)
                    nuisance(:,ii) = smooth(nuisance(:,ii));
                    dtnuisance(:,ii) = smooth(dtnuisance(:,ii));
                    sqnuisance(:,ii) = smooth(sqnuisance(:,ii));
                end
            end
            motionDer = [nuisance dtnuisance sqnuisance];
            %can
            if useGlobal
                if useNuisance
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
            else
                if useNuisance
                    if useNuisGas
                        np = []; np = [motionDer L' nuis_res'];
                    else
                        np = []; np = [motionDer L'];
                    end
                else
                    if useNuisGas
                        np = []; np = nuiss_probe';
                    else
                        np = [];
                    end
                end
            end
            
            %remove signals explained by nuissance regressors
            [cleanBOLD] = scrubBOLD(rBOLD, mWBmask, nprobe', np, opts.figdir);
            
            if denoise
                if license('test', 'wavelet_toolbox')
                    [ cleanBOLD ~ ]  = wavDenRO1(cleanBOLD, mWBmask);
                    saveas(gcf,[opts.figdir,'dBOLD.fig']);
                else
                    %perform some temporal smoothing of data
                end
            end
            
            if smthdata
                cleanBOLD( cleanBOLD == 0) = NaN;
                classname = class( cleanBOLD);
                mWBmask = eval([classname,'(mWBmask)']);    
                
                cleanBOLD = smthData( cleanBOLD, mWBmask, opts);
            end
            
            cleanBOLD = normTimeseriesMM(cleanBOLD, mWBmask, [5 30]);
            
            if save_cleaned; saveImageData(cleanBOLD, BOLDheader, opts.resultsdir, 'cleanBOLD.nii.gz', 64); end
            
            %setup lag analysis directories
            opts.corrlagdir = [opts.resultsdir,'CORRlag\']; mkdir(opts.corrlagdir);
            opts.glmlagdir = [opts.resultsdir,'GLMlag\']; mkdir(opts.glmlagdir);
            opts.loadprobe = 0;
            BOLDsig_analysis(GMmask.*mWBmask,mWBmask,cleanBOLD,nprobe',nuiss_probe', opts)
            clear voxel_ts error_data np_coef nuisance dtnuisance cleanBOLD nnuiss nnnuiss nnuis_TS nnnuis_TS probe nuis_probe
            save([opts.resultsdir,'processing_options.mat'], 'opts');
        end
    end
end

%% Generate CVR, lag, corrected CVR maps (correlation & GLM based)

loadprobe = 1; %use BOLD regressor generated above instead of endtidal trace
smoothTempHRF =1;
HRF_mapsRO1(TR,headers,GMmask.*mWBmask,mWBmask,BOLD_ts,probe(1,x_ind(1):x_ind(2))',smoothTempHRF,prewhite,interp_factor,dim,savedir1,nuiss_probe',loadprobe)
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
