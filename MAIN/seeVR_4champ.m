%seeVR code writen by Alex Bhogal (a.bhogal@umcutrecht.nl). Do not share
%code without express permission from A. Bhogal

clear all
close all
clear global all

%4 THE Champ
%set your own paths pointing to BOLD data and breathing files
subj = 'bh04'
task = 'exp'
dir = ['C:\Users\abhogal\Documents\DATA\MRI_data\Champagne\',subj,'\mri\BOLD_',task,'\'];
seqpath  = ['C:\Users\abhogal\Documents\DATA\MRI_data\Champagne\',subj,'\ramr\bh.exp.x4.bold.mb.',task,'\'];
%seqpath  = ['C:\Users\abhogal\Documents\DATA\MRI_data\Champagne\',subj,'\ramr\hc.step.bold.mb\'];

cd(dir)
%load motion/distortion corrected BOLD data
filename = ls('*BOLD_applytopup.nii.gz*')
[BOLD,INFO,BOLDheader] = loadImageData(dir, filename);
%load GM segmentation
file = ls('*seg_1*');
[GMmask,INFOmask,HEADERmask] = loadImageData(dir, file);
GMmask(isinf(GMmask)) = 0; GMmask = double (GMmask);
%load WM segmentation
file = ls('*seg_2*');
[WMmask,~,~] = loadImageData(dir, file);
WMmask(isinf(WMmask)) = 0; WMmask = double (WMmask);
%load CSF segmentation
file = ls('*seg_0*');
[CSFmask,~,~] = loadImageData(dir, file);
CSFmask(isinf(CSFmask)) = 0;  CSFmask = double (CSFmask);
file = ls('brain_mask*');
[WBmask,~,~] = loadImageData(dir, file);
WBmask(isinf(WBmask)) = 0; WBmask = double (WBmask);


%%%miscellaneous options for directories, data, loading
global opts
opts.TR = BOLDheader.dime.pixdim(5);
opts.voxelsize = BOLDheader.dime.pixdim(2:4);
opts.dyn = size(BOLD,4); xdata = [opts.TR:opts.TR:opts.TR*opts.dyn];
opts.savedir = dir;
opts.figdir = [dir,'figures\']; mkdir(opts.figdir);
opts.headers = struct();
opts.headers.ts = BOLDheader;
opts.headers.map = opts.headers.ts;
opts.headers.map.dime.dim(5) = 1;  
opts.headers.mask = HEADERmask;
%% load breathing traces
%%%options for outlier removal
opts.seqpath = uigetdir(seqpath) %select main directory for session
opts.extra = 30; %defines number of extra samples taken before and after the start and end of the events
opts.remOUT = 1; %removes ourliers in general (both positive and negative). Doesnt always work so take care.
opts.remOUTbh = 1; %removes only outliers lower than baseline. This is specifically for BH as we dont expect low PetCO2. Also set 'opts.remOUT = 1' to work
%you can modify the events files so the load functions finds the start
%and end points of the sequence. For now I turned that feature off for
%the BH
opts.evenstart = 'SequenceStart';
opts.eventend = 'SequenceEnd';
[corrvec_CO2,corrvec_O2] = loadRAMRgen4(opts); %loads RAMR gen4 gas traces
CO2_corr = corrvec_CO2;
O2_corr = corrvec_O2;
x_indyi(1) = 1; x_indyi(2)= length(CO2_corr);

%% Align gas traces with BOLD data
cd(dir)
%generate timeseries for correlations
% convolve CO2 trace with HRF
opts.HRF = 1;
opts.onset = 1;
opts.disp = 20
opts.under = 1;
opts.plot = 1;
[~,HRF_CO2_corr] = convHRF(CO2_corr, opts);

SE = strel('diamond',1);
TS = meanTimeseries(BOLD,GMmask);
%GUI to correct alignment
trAlign(CO2_corr, O2_corr, TS, opts);
CO2trace = probe1;
O2trace = probe2;

%% generate global signal and perform initial data 'scrubbing'
cd(dir)
%%%preprocessing options

opts.motioncorr = 0.7; %motion parameters with higher correlation than threshold will not be included in the regression
opts.legOrder = [0]; %order of Legendre polynomials to include in regression (typically up to 4th order however can cause problems)
opts.plot = 0;
opts.HRF = 1;

opts.uni = 0;% if this is set to 1 negative correlations will be ignored
opts.loadprobe = 0; %saves time for creating regressor if one with the correct length already exists

%%%refinement optimized regressor options
opts.corr_thresh = 0.5; %threshold by which to accept correlated voxels during the refinement of the BOLD regressor
opts.save_cleaned = 0; %set to 1 to save BOLD data with motion/legendre regressed out



%% pre-process data and run preliminary analysis
taskbh = 0;
for bhidx = [240 370 530 size(BOLD,4)] 
    taskbh = taskbh+1
opts.priority = [1] %priority = aa %1 for CO2, 0 for O2
opts.denoise =  [0] %perform wavelet denoising in the temporal dimension
for useGlobal = [1 0] %use global signal regressor during data-cleanup


            opts.smthdata = 1; %turn on spatial smoothing
            
            cd(dir)
            if opts.priority
                probe = CO2trace;
                opts.useNuisGas = 0; %use alternative gas as nuissance (most important when looking at O2 as primary stimulus
                nuiss_probe = O2trace;
                opts.prefix = [task,int2str(taskbh)]
            else
                probe = O2trace;
                opts.useNuisGas = 1;
                nuiss_probe = CO2trace;
                opts.prefix = [task,int2str(taskbh)]
                
            end
            
            %update figure and results directories
            opts.resultsdir = [opts.savedir,opts.prefix,'_GS',int2str(useGlobal),'_WD',int2str(opts.denoise),'_CO2pr',int2str(opts.priority),'\'];
            mkdir(opts.resultsdir)
            opts.figdir = [opts.resultsdir,'figures\']; mkdir(opts.figdir);
            
            %reduce BOLD series to conserve memory for processing later on:
            %3rd argument can be start/end indices; if blank, user selects
            [idx, rBOLD] = chopTimeseries(BOLD,WBmask, [30 bhidx]);
            %update options variables based on new epoch
            nuiss_probe = rescale(nuiss_probe(1,idx(1):idx(2)));
            nprobe = []; nprobe = probe(1,idx(1):idx(2));
            CO2_probe = CO2trace(1,idx(1):idx(2));
            O2_probe = O2trace(1,idx(1):idx(2));
            
            if opts.HRF
            %parameter vetors can be adjusted based on the size/TR of the
            %data
            %Convolve with double gamma tuned for onset
            opts.onset = 1:1*30;
            opts.disp = 0.4;
            opts.under = 1;
            [HRF,onHRFprobe] = convHRF(nprobe,opts); %Convolve with double gamma tuned for onset
            
            normprobe = onHRFprobe(1,:);
            onHRFprobe(1,:) = [];
            cshift = -round(5/opts.TR);
            for ii=1:size(onHRFprobe,1); onHRFprobe(ii,:) = circshift(onHRFprobe(ii,:),cshift); onHRFprobe(ii,end+cshift:end) = onHRFprobe(ii,end+cshift); end
            figure(30); subplot(1,2,1); plot(onHRFprobe'); hold on; plot(normprobe', 'k', 'LineWidth', 2); title('Shifted onset, limited dispersion')
            
            %Convolve with double gamma tuned for dispersion
            opts.onset = 1;
            opts.disp = 1:1*20
            opts.under = 1;
            [HRF,dspHRFprobe] = convHRF(nprobe, opts);
            
            %removes probe (since we already have a probe vector)
            dspHRFprobe(1,:) = [];
            figure(30); subplot(1,2,2); plot(dspHRFprobe'); hold on; plot(normprobe', 'k', 'LineWidth', 2); title('Dispersion')
            
            HRFprobe = [ onHRFprobe' dspHRFprobe' ]
            else
                HRFprobe = [];
            end
            %load motion parameters (based on FSL, may need to be adjusted
            %for other MP files)
            mpfilename = ls('*mcf.par');
            nuisance = load(mpfilename); nuisance = nuisance(idx(1):idx(2),:);
            dtnuisance =  gradient(nuisance);
            sqnuisance = nuisance.*nuisance;
            motionDer =[nuisance dtnuisance sqnuisance];
            
            %initialize Legendre regressors
            if opts.legOrder == 0; L = []; else
                L = []; for ii=opts.legOrder; L(ii,:) = rescale(LegendreN(ii,opts.xdata)); end
            end
            %setup nuisance regressors to be used with HRFprobe to
            %determine global signal residual (for O2 priority, CO2 should
            %always be regressed out)
            if opts.useNuisGas
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
            [~,nuis_res,resBOLD] = genGS(rBOLD, WBmask, np,HRFprobe, opts);
            
            clear resBOLD
            
            % Final cleanup of BOLD data using specified nuissance regressors
            %remove nuisance correlating with probe
            autoCorr = abs(corr(meanTimeseries(rBOLD,mWBmask)',motionDer)); %remove
            autoCorr(autoCorr < opts.motioncorr) = 0; autoCorr(autoCorr > 0) = 1; %removes anything with more than weak correlation
            %remove highly correlated nuisance regressors to preserve signal response
            index = ([1:1:size(motionDer,2)]).*autoCorr; index(index == 0) = [];
            motionDer(:,index) = [];
            
            %setup nuisance regressors - explore the effect of adding GSR
            if useGlobal
                if opts.useNuisGas
                    np = []; np = [motionDer L' nuis_res' nuiss_probe'];
                else
                    np = []; np = [motionDer L'  nuis_res'];
                end
            else
                if opts.useNuisGas
                    np = []; np = [motionDer L' nuiss_probe'];
                else
                    np = []; np = [motionDer L' ];
                    
                end
            end
            [np,~]=licols(np); %remove linearly dependent components
            
            %remove signals explained by nuissance regressors
            [scleanBOLD] = scrubData(rBOLD, WBmask, np, nprobe, opts.figdir);
            if opts.denoise; [cleanBOLD] = denoiseData(cleanBOLD,WBmask,opts); end
            %spatial smoothing
            opts.spatialdim = 3; opts.filtWidth = 7; opts.FWHM = 3; 
            cleanBOLD = smthData( scleanBOLD, WBmask, opts);
                          
            %normalize processed timeseries to baseline (if no index
            %supplied in third argument, user may choose begin/end baseline
            cleanBOLD = normTimeseries(cleanBOLD, WBmask, [5 25]);
            %save cleaned timesieres
            if opts.save_cleaned; saveImageData(cleanBOLD, BOLDheader, opts.resultsdir, 'cleanBOLD.nii.gz', 64); end
            
            %setup lag analysis directories
            %2 options to consider: 1) supply 'clean' data for correlation
            %analysis (and shifter GLM) OR 2) supply original data +
            %nuissance regressors for GLM analysis.
            %%%%%% 1)
            opts.lowerlagthresh = -5
            opts.upperlagthresh = 5
            opts.interp_factor = 4;
            opts.corr_thresh = 0.5;
            opts.glm_model = 0;
            opts.uni = 1;
            opts.lowlag = -5
            opts.highlag = 20;
            SE = strel('diamond',1);
            lagCVR(imerode(GMmask,SE).*mWBmask,mWBmask,cleanBOLD,nprobe,nuiss_probe, opts); %scrubbed data

            %save options for reference
            opts.bhidx = bhidx;
            save([opts.resultsdir,'processing_options.mat'], 'opts');

end
end