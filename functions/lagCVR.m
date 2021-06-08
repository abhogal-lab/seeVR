 %% PCA, RAPIDTIDE/RIPTIDE by Allen Champagne - CNS, Queen's University, ON on July 16th 2017
%% Expanded by Alex Bhogal - UMC Utrecht

% All copyrights are reserved to the authors. Please do not share this code
% without the permission of the authors.

% This script estimates the lag and (correctd)CVR maps and provides statistical outputs
% using a modified version of the RIPTIDE method (Donahue 2016). This
% Outputs can be created using a correlation or GLM based approach


function [newprobe] = lagCVR(GMmask,mask,BOLD_ts,probe,nuisance,opts)
global opts; 

if iscolumn(probe); else; probe = probe'; end
test1 = nuisance(1,:); test2 = nuisance(:,1);
if length(test1) > length(test2); nuisance = nuisance'; end; clear test1 test2

%setup default parameters
if isfield(opts,'prewhite'); else; opts.prewhite = 0; end                 %zero mean and unit variance of data. normally leads to bad results
if isfield(opts,'interp_factor'); else; opts.interp_factor = 4; end       %factor by which to temporally interpolate data. Better for picking up lags between TR
if isfield(opts,'load_probe'); else; opts.load_probe = 0; end             %saves time for creating regressor if one with the correct length already exists
if isfield(opts,'save_rts'); else; opts.save_rts = 0; end                 %save correlation timeseries - can be used to visualize lags
%setup main default parameters
if isfield(opts,'rescale_probe'); else; opts.rescale_probe = 1; end       %rescaling may be helful for initial refinement run
if isfield(opts,'trace_corr'); else; opts.trace_corr = 1; end             %perform additional correlation with gas trace (on top of optimized regressor)
if isfield(opts,'refine_regressor'); else; opts.refine_regressor = 1; end %refine BOLD regressor. If '0' a straight correlation will be done with probe
if isfield(opts,'pca_analysis'); else; opts.pca_analysis = 1; end         %PCA analysis to optimize BOLD regressor
if isfield(opts,'corr_model'); else; opts.corr_model = 1; end             %perform correlation based analysis
if isfield(opts,'cvr_maps'); else; opts.cvr_maps = 1; end                 %generate CVR maps based on regression of BOLD signal with gas probe
if isfield(opts,'eff_probe'); else; opts.eff_probe = 1; end               %uses effective probe to calculate probe response
if isfield(opts,'glm_model'); else; opts.glm_model = 0; end               %perform GLM based lag regression. This works well for high temporal resolution data
if isfield(opts,'uni'); else; opts.uni = 0; end                           %if this is set to 1 negative correlations will be ignored
if isfield(opts,'norm_regr'); else; opts.norm_regr = 0; end               %normalize the regressor for correlation analysis between 0-1
%important parameters for results
if isfield(opts,'corr_thresh'); else; opts.corr_thresh = 0.7; end         %threshold by which to accept correlated voxels during the refinement of the BOLD regressor
%refinement
if isfield(opts,'lowerlagthresh'); else; opts.lowerlagthresh = -2; end    %lower threshold for refinement (generaly -3 to 0)
if isfield(opts,'upperlagthresh'); else; opts.upperlagthresh = 2; end     %upper threshold for refinement (generaly 0 to +3)
%correlation
if isfield(opts,'lowlag'); else; opts.lowlag = -5; end                    %lower threshold for correlation (generaly -3 to 0)
if isfield(opts,'highlag'); else; opts.highlag = 90; end                  %upper threshold for correlation (in healthy up to 20-40, in disease 60-90)

%account for interpolation factor
opts.adjlowlag = round((opts.lowlag/opts.TR)*opts.interp_factor); %setup lower lag limit; negative for misalignment and noisy correlation
opts.adjhighlag = round((opts.highlag/opts.TR)*opts.interp_factor); %setups upper lag limit; allow for long lags associated with pathology
opts.adjlowerlag = round(opts.lowerlagthresh/opts.TR*opts.interp_factor);
opts.adjupperlag = round(opts.upperlagthresh/opts.TR*opts.interp_factor);

%check
if opts.load_probe == 0 && opts.refine_regressor == 0 && opts.trace_corr == 0
    disp('check options; stopping analysis')
    opts.glm_model = 0;
    opts.corr_model = 0;
    opts.cvr_maps = 0;
end

%setup save directories
if opts.corr_model; opts.corrlagdir = [opts.resultsdir,'corrLAG\']; mkdir(opts.corrlagdir); end
if opts.corr_model && opts.cvr_maps; opts.corrCVRdir = [opts.corrlagdir,'CVR\']; mkdir(opts.corrCVRdir); end
if opts.corr_model == 0; opts.cvr_maps = 0; end
if opts.glm_model 
    opts.glmlagdir = [opts.resultsdir,'glmLAG\']; mkdir(opts.glmlagdir);
    opts.glmCVRdir = [opts.glmlagdir,'CVR\']; mkdir(opts.glmCVRdir);
end
%optinal image outputs
if isfield(opts, 'robustTstat'); else; opts.robustTstat = 1; end
if isfield(opts, 'robustR'); else; opts.robustR = 0; end

cd(opts.resultsdir);
datatype = 64;
[xx, yy, zz, dyn] = size(BOLD_ts);
orig_regr = probe;

t = cputime;

%% grab coordinates
% WB coordinates
[orig_voxel_ts, coordinates] = grabTimeseries(BOLD_ts, mask);
% GM coordinates
[gm_voxel_ts, gmcoordinates] = grabTimeseries(BOLD_ts, GMmask);

if opts.prewhite
    gm_voxel_ts = gm_voxel_ts'; parfor ii=1:size(gmcoordinates,1); [gm_voxel_ts(:,ii), ~, ~] = prewhiten(gm_voxel_ts(:,ii)); end; gm_voxel_ts = gm_voxel_ts';
    pw_voxel_ts = orig_voxel_ts'; parfor ii=1:size(coordinates,1); [pw_voxel_ts(:,ii), ~, ~] = prewhiten(pw_voxel_ts(:,ii)); end; pw_voxel_ts = pw_voxel_ts';
    [probe, ~, ~] = prewhiten(probe);
end

%interpolate variables
probe = interp(probe,opts.interp_factor); %interpolate data by a factor of the number of slices for more accurate timing
orig_regr = interp(orig_regr,opts.interp_factor); %save input regressor to generate corrected CVR maps
if opts.glm_model; for ii=1:size(nuisance,2); np(:,ii) = interp(nuisance(:,ii),opts.interp_factor); end; end

gm_voxel_ts_nonan=zeros([length(gmcoordinates),opts.interp_factor*dyn]);
parfor ii = 1:length(gmcoordinates)
    gm_voxel_ts_nonan(ii,:) = interp(gm_voxel_ts(ii,:),opts.interp_factor);
end

%% to save time, check if there is an existing probe
if opts.refine_regressor
    probename = 'final_probe.mat'
    if exist(probename) && opts.load_probe
        newprobe = load(probename); newprobe = newprobe.newprobe;
        disp('found existing probe; skipping generation of BOLD regressor')
    else
        
        %% RipTIDE based on Donahue/Frederick 2016
        
        disp('Creating optimized regressor')
        %grab GM timeseries to generate probe
        
        keep_probes = [];
        %keeps looping until RSME between probe is very small = therefore, reach
        %convergence!
        roundprobe=0;
        stop=0;
        while stop==0
            roundprobe=roundprobe+1;
            
            if roundprobe == 1 %rescaling probe may be helpful for initial iteration
                if opts.rescale_probe
                newprobe = rescale(probe);
                else
                newprobe = probe;
                end
            end
            
            disp(['Correlating TS voxelwise to probe ',num2str(roundprobe)]);disp(' . ');disp(' . ');disp(' . ');
            % voxelwise corrrelation with probe
            
            [~,lags] = xcorr(newprobe,gm_voxel_ts_nonan(1,:),'coeff');  %%%%%% !!!!!!!!!!
            %matrix cross correlation with probe
            tic
            a=gm_voxel_ts_nonan;
            a2=mat2cell(a,ones(1,size(a,1)),size(a,2)); %turn into cells so treat all rows independently
            
            % the order here matters because of the SIGN of the lag
            % this way we slide the probe to the ts
            b2=cellfun(@(a2) xcorr(a2,newprobe,'coeff'),a2,'uniformoutput',false);
            corr_probe=cell2mat(b2); %reformat into correlation matrix
            corr_probe = corr_probe' ;
            
            %remove low and high lags ignoring unreasonable and long lag times
            idx = find(lags<=opts.adjlowlag | lags>=opts.adjhighlag);
            lags(idx)=[];
            corr_probe(idx,:)=[];      clear idx
            
            % estimate peak
            % delay and correlation over specified time lag range with fit parameters
            maxcorr = max(corr_probe,[],1);
            tmp = repmat(maxcorr,length(corr_probe(:,1)),1);
            yes = zeros(size(corr_probe));
            yes(find(tmp==corr_probe))=1; clear tmp
            yes(yes==0)=NaN;
            tmp = repmat(lags',1,length(corr_probe(1,:)));
            ts_lag = nanmean(tmp.*yes); %clean things up
            
            checkcorr=maxcorr;
            %creating a filter to remove the timeseries before the PCA
            checkcorr(maxcorr>=opts.corr_thresh)=0;  % if good, put 0 because will be kept
            checkcorr(maxcorr<opts.corr_thresh)=1;  % to be removed
            
            checklag=ts_lag;
            checklag(ts_lag>=opts.adjlowerlag & ts_lag<=opts.adjupperlag)=0;
            checklag(ts_lag<opts.adjlowerlag | ts_lag>opts.adjupperlag)=1;
            
            resultant_filter = zeros(size(checklag));
            resultant_filter(checklag==1 | checkcorr==1)=1;
            remove = find(resultant_filter==1);
            
            % removing timeseries and lags outside of range
            gm_voxel_ts_nonan_filt=gm_voxel_ts_nonan';
            gm_voxel_ts_nonan_filt(:,remove)=[];
            ts_lag_filt=ts_lag;
            ts_lag_filt(remove)=[];
            
            % re-aligning all timeseries based on lag for PCA
            xfit = 1:opts.interp_factor*dyn;
            
            new_shifted_ts = zeros([length(ts_lag_filt), opts.interp_factor*dyn]);
            parfor hh = 1:length(ts_lag_filt)
                new_shifted_ts(hh,:) = interp1(xfit,gm_voxel_ts_nonan_filt(:,hh),xfit+ts_lag_filt(hh),'Linear',gm_voxel_ts_nonan(hh,1));
            end
            
            disp('finished correlation, estimating new probe')
            %license('test','Statistics_toolbox')
            if opts.pca_analysis
                
                %%% PCA on significant shifted timeseries to get single timecourse that explains
                %the highest shared variance
                % score = how much pattern from each component is in the voxel
                % coeff = actual component pattern
                
                [coeff,score,latent,~,explained] = pca(new_shifted_ts,'VariableWeights','variance','Row','Complete'); % complete option ignores NaN
                
                %loop until have enough components to explain 85% of the
                %variability
                total_components = 3;jj=total_components-1;  % minimum 3 components
                while jj < length(explained)
                    jj=jj+1;
                    checkk = sum(explained(1:jj));
                    if checkk >= 85 % if first 3 component covers more than 85, continue
                        continue
                    else
                        total_components = total_components+1; % see explained, usually explains ~ 75-85% with first 10 components
                    end
                end
                
                try
                pca_components = score(:,1:total_components)*coeff(:,1:total_components)';
                catch 
                disp('The correlation parameter for the optimized regressor is set to high')
                disp('i.e. there are not sufficient highly correlated voxels: set opts.corrthresh to a lower value')
                error('exiting function')
                return
                end
                pca_components = pca_components(~any(isnan(pca_components),2),:); %removing rows with NaN
                
                % components can be flipped so doing a pearson correlation with the
                % probe gas trace; if flipped, will have a - correlation so
                % just multiply by -1
                outputpca = nan(size(pca_components));
                clear pearsR
                pearsR = zeros([1 length(pca_components(:,1))]);
                tic
                parfor dd = 1:length(pca_components(:,1))
                    pearsR(dd) = corr(newprobe,pca_components(dd,:)');
                    
                    if pearsR(dd) > 0
                        outputpca(dd,:) = 1 * pca_components(dd,:);
                    else
                        outputpca(dd,:) = -1 * pca_components(dd,:);
                    end
                    
                end
                toc
                % weigthed average of refined regressor
                newprobe=[];
                newprobe = nanmean(outputpca,1)';
                
            else
                newprobe=[];
                newprobe = (nanmean(new_shifted_ts,1))';
            end
            
            % averaging all targeted components for new probe
            keep_probes(roundprobe,:) = newprobe;
            clear max_corr_voxel   pca_components  outputpca   new_shifted_ts pearsR
            
            rmse=[];
            if roundprobe==1
                rmse = immse(probe(:),newprobe(:))
            else
                rmse = immse(newprobe',(keep_probes(roundprobe-1,:)))
            end
            
            if rmse>0 & rmse<0.005
                stop = 1;
            elseif roundprobe==10 %maximum number of iterations to find probe
                stop = 1;
            end
            
        end
        figure(15);
        plot(keep_probes','LineWidth',2);title('All probes');
        saveas(gcf,[opts.resultsdir,'all_probes.fig']);
        %save the final probe
        save([opts.resultsdir,'final_probe.mat'], 'newprobe');
        
        disp('Finished creating optimized regressor')
        
        clear a GMmask yes score coeff gm_voxel_ts gm_voxel_ts_nonan_filt
    end
else
    newprobe = probe;
    opts.trace_corr = 0;
    [~,lags] = xcorr(newprobe,gm_voxel_ts_nonan(1,:),'coeff');
    idx = find(lags<=opts.adjlowlag | lags>=opts.adjhighlag);
    lags(idx)=[];
end
%% for final round, take the refined regressor and do final
% correlations between it and the original signal to estimate the peak time delay,
% correlation for the whole brain

%%% grabbing WB timeseries
tic
wb_voxel_ts=zeros([length(coordinates) opts.interp_factor*dyn]);
clean_voxel_ts=zeros([length(coordinates) opts.interp_factor*dyn]);
%For lags you want to use wb_voxel_ts (accounting for prewhitening. For CVR
%maps always regress original regressor with clean timeseries
ip = opts.interp_factor;
if opts.prewhite
    parfor ii = 1:length(coordinates)
        clean_voxel_ts(ii,:) = interp(orig_voxel_ts(ii,:),ip);
        wb_voxel_ts(ii,:) = interp(pw_voxel_ts(ii,:),ip);
    end
else
    parfor ii = 1:length(coordinates)
        clean_voxel_ts(ii,:) = interp(orig_voxel_ts(ii,:),ip);
    end
    wb_voxel_ts = clean_voxel_ts;
end; clear ip


if opts.corr_model
    
    if opts.trace_corr; qq = [1 2]; else; qq = 1; end %correlate probe y/n
    for pp = qq
        switch pp
            case 1
                %setup regression  using optimized regressor
                regr = newprobe;
                disp('performing correlation based lag analysis using optimized regressor')

            case 2
                %setup regression using input probe
                disp('performing correlation based lag analysis using input probe')
                regr = probe;
        end
        
        % find lag between probe and timeseries
        corr_probe=[];idx_matrix=[]; a2=[];
        %matrix cross correlation with probe
        a2=mat2cell(wb_voxel_ts,ones(1,size(wb_voxel_ts,1)),size(wb_voxel_ts,2)); %turn into cells so treat all rows independently
        %THE ORDER HERE MATTERS!!!
        b2=cellfun(@(a2) xcorr(a2,regr,opts.highlag,'normalized'),a2,'uniformoutput',false);
        corr_probe=cell2mat(b2); %reformat into correlation matrix
        corr_probe = corr_probe';
        %use highlag to restrict the number of correlations that need to be done
        [~,lags] = xcorr(gm_voxel_ts_nonan(1,:),regr,opts.highlag,'normalized');  %%%%%% !!!!!!!!!!
        
        %save correlations overtime
        corr_ts = zeros([xx*yy*zz size(corr_probe,1)]);
        corr_ts(coordinates,:) = corr_probe';
        corr_ts(isnan(corr_ts)) = 0;
        corr_ts = reshape(corr_ts, [xx yy zz size(corr_probe,1)]);
        % trim correlation probe
        idx = find(lags<=opts.adjlowlag | lags>=opts.adjhighlag);
        %remove lag values lower that what is possible
        lags(idx)=[];
        corr_probe(idx,:)=[];      clear idx
        
        %create lag maps
        rvec = zeros([1, length(coordinates)]);
        index = zeros([1 size(wb_voxel_ts,1)]); tmpmaxcorr = zeros([1 size(wb_voxel_ts,1)]);
        
        if opts.uni
            [~, index] = max(corr_probe); %negative correlations will not weigh more than positive ones
        else
            [~, index] = max(abs(corr_probe)); %Absolute max to take care of negative correlation
        end
        %fill correlation map
        parfor ii = 1:length(coordinates)
            rvec(1,ii) = corr_probe(index(1,ii),ii);
        end
        
        r_map = zeros([xx*yy*zz 1]); lag_map = zeros([1 xx*yy*zz]);
        r_map(coordinates) = rvec;
        r_map = reshape(r_map,[xx yy zz]);
        %fill lag map
        lag_map(1,coordinates) = lags(1,index);
        tmpLag = opts.TR*(lag_map/opts.interp_factor); tmpLag(1,coordinates) = tmpLag(1,coordinates) + abs(min(tmpLag(:)));
        lag_map = reshape(lag_map,[xx yy zz]);
        tmpLag = reshape(tmpLag,[xx yy zz]);
        %save lag and r maps
        
        switch pp
            case 1
                saveImageData(mask.*tmpLag, opts.headers.map, opts.corrlagdir, 'lag_map.nii.gz', datatype)
                saveImageData(mask.*lag_map, opts.headers.map, opts.corrlagdir, 'uncorr_lag_map.nii.gz', datatype)
                saveImageData(mask.*r_map, opts.headers.map, opts.corrlagdir, 'r_map.nii.gz', datatype)
 
                if opts.save_rts
                    disp('saving correlation timeseries based on optimized regressor')
                    saveImageData(corr_ts, headers.ts, opts.corrlagdir,  'r_ts.nii.gz', datatype)
                clear Trcorr_ts
                end
                LAG(1,:,:,:) = mask.*mask.*tmpLag;
                RL(1,:,:,:) = mask.*r_map;
                clear tmpLag lag_map r_map
            case 2
                saveImageData(mask.*tmpLag, opts.headers.map, opts.corrlagdir,  'lag_map_probe.nii.gz', datatype);
                saveImageData(mask.*lag_map, opts.headers.map,opts.corrlagdir, 'uncorr_lag_map_probe.nii.gz', datatype);
                saveImageData(mask.*r_map, opts.headers.map, opts.corrlagdir, 'r_map_probe.nii.gz', datatype);
                if opts.save_rts
                    disp('saving correlation timeseries based on input probe')
                    saveImageData(Trcorr_ts, opts.headers.ts, opts.corrlagdir, 'r_ts_probe.nii.gz', datatype);
                clear Trcorr_ts
                end
                LAG(2,:,:,:) = mask.*mask.*tmpLag;
                RL(2,:,:,:) = mask.*r_map;
                clear tmpLag lag_map r_map
        end
        clear Trcorr_ts orig_voxel_ts
    end
    
    disp(['Lag, regressor and r_maps were created in: ',int2str((cputime-t)/60),' minutes'])
    
    clear gm_voxel_ts_nonan m corr_ts regr
end


%% if correlating both input regressor and optimized probe
if opts.trace_corr
    disp('Calculating robust lag map')
    [MR,IR] = max(abs(RL),[],1,'omitnan'); IR = squeeze(IR);
    LAG = reshape(LAG,2,xx*yy*zz);
    IR = reshape(IR,1,numel(IR));
    robustIR = zeros(size(IR));
    for ii=1:length(IR)
        tmp = IR(1,ii);
        robustIR(1,ii) = LAG(tmp,ii);
    end
    robustIR = reshape( robustIR,size(mask));
    %save image
    saveImageData(mask.*robustIR, opts.headers.map, opts.corrlagdir,'robustLAG_r.nii.gz', datatype);
    clear LAG RL
else
    clear LAG RL
end


    %% Generate lag adjusted CVR maps

if opts.cvr_maps
    %% normalize probe and refined regressor, determine effective CO2 trace
    
    if opts.eff_probe; qq = [1 2]; else; qq = 1; end %regress with probe y/n
    for pp = qq
        switch pp
            case 1
                %setup regression  using optimized regressor
                disp('Generating base maps using probe trace')
                regr = probe;
            case 2
                %setup regression using end-tidal CO2
                rs_newprobe = rescale(newprobe, 0.001,1 );
                coef = glmfit(rs_newprobe, probe);
                eff_probe = coef(2,1)*rs_newprobe + coef(1,1);
                
                figure; subplot(3,1,1); plot(probe,'b'); title('initial probe probe'); ylabel('mmHg')
                subplot(3,1,2); plot(newprobe); title('BOLD regressor'); ylabel('% BOLD change')
                subplot(3,1,3); plot(probe,'b'); hold on; plot(eff_probe, 'r'); title('effective probe'); ylabel('mmHg')
                saveas(gcf,[opts.figdir,'effective_probes.fig']);
                disp('Generating base maps using effective probe trace')
                regr = eff_probe;
        end
        
        
        bCVR = zeros([1 xx*yy*zz]); bR2 = zeros([1 xx*yy*zz]); bSSE = zeros([1 xx*yy*zz]); bTstat = zeros([1 xx*yy*zz]);
        
        %create base CVR map
        A = regr;
        C = [ones([length(A) 1]) A]; clear A
        regr_coef = C\clean_voxel_ts'; %perform least squares linear regression
        
        bCVR(1,coordinates) = regr_coef(2,:); %extract slope
        bCVR(bCVR > 20) = 0; bCVR(bCVR < -20) = 0; %cleanup base CVR map
        bCVR = reshape(bCVR, [xx yy zz]);
        %original observations
        X = wb_voxel_ts;
        %fitted values
        Y = regr_coef(2,:).*repmat(regr,[1 size(X,1)]) + ones([1 length(regr_coef)]).*regr_coef(1,:);
        
        SSE = zeros([1 size(X,1)]); SST = zeros([1 size(X,1)]); STDEVr = zeros([1 size(X,1)]);
        parfor ii=1:size(X,1)
            STDEVr(1,ii) = nanstd(X(ii,:)-Y(:,ii)'); %standard deviation of residuals (t = beta/STDEVr)
            SSE(1,ii) = (norm(X(ii,:) - Y(:,ii)'))^2;
            SST(1,ii) = (norm(X(ii,:)-mean(X(ii,:))))^2;
        end
        
        bT = regr_coef(2,:)./STDEVr(1,:);
        R2 = 1 - SSE./SST;
        bR2(1, coordinates) = R2; bR2 = reshape(bR2, [xx yy zz]);
        bSSE(1, coordinates) = SSE; bSSE = reshape(bSSE, [xx yy zz]);
        bTstat(1, coordinates) = bT; bTstat = reshape(bTstat, [xx yy zz]);
        switch pp
            case 1 %entdidal
                %save base CVR map
                saveImageData(mask.*bCVR, opts.headers.map, opts.corrCVRdir, 'bCVR_map.nii.gz', datatype);
                %save base stats data
                saveImageData(mask.*bR2, opts.headers.map, opts.corrCVRdir,'bR2_map.nii.gz', datatype);
                saveImageData(mask.*bSSE, opts.headers.map, opts.corrCVRdir,'bSSE_map.nii.gz', datatype);
                saveImageData(mask.*bTstat, opts.headers.map, opts.corrCVRdir,'bTstat_map.nii.gz', datatype);
                disp('Base CVR, r2, Tstat and SSE maps were created using entidal regressor')
                CVR(1,:,:,:) = mask.*bCVR;
                TSTAT(1,:,:,:) = mask.*bTstat;
                RC(1,:,:,:) = mask.*bR2;
                clear bCVR bR2 bT bTstat bSSE SSE SST R2 X Y r regr_coef A C s
            case 2 %effective probe
                %save base CVR map
                saveImageData(mask.*bCVR, opts.headers.map, opts.corrCVRdir, 'bCVR_eff_map.nii.gz', datatype);
                %save base stats data
                saveImageData(mask.*bR2, opts.headers.map, opts.corrCVRdir,'bR2_eff_map.nii.gz', datatype);
                saveImageData(mask.*bSSE, opts.headers.map, opts.corrCVRdir,'bSSE_eff_map.nii.gz', datatype);
                saveImageData(mask.*bTstat, opts.headers.map, opts.corrCVRdir,'bTstat_eff_map.nii.gz', datatype);
                disp('Base CVR, r2, Tstat and SSE maps were created using effective entidal regressor')
                CVR(2,:,:,:) = mask.*bCVR;
                TSTAT(2,:,:,:) = mask.*bTstat;
                RC(2,:,:,:) = mask.*bR2;
                clear bCVR bR2 bT bTstat bSSE SSE SST R2 X Y r regr_coef A C s
        end
        
        
        %% generate lag-corrected CVR maps
        
        disp('Generating lag-adjusted maps')
        
        
        cCVR = zeros([1 xx*yy*zz]); cR2 = zeros([1 xx*yy*zz]); cSSE = zeros([1 xx*yy*zz]); cTstat = zeros([1 xx*yy*zz]);
        
        shifted_regr = NaN([length(coordinates) dyn*opts.interp_factor]);
        
        for ii = 1:length(coordinates)
            corr_regr = circshift(regr',index(ii));
            if index > 0
                corr_regr(1,1:index(ii)) = NaN;
            elseif index < 0
                corr_regr(1,end-index:end) = NaN;
            else
                %do nothing because index is zero = zero shift
            end
            shifted_regr(ii,:) = corr_regr;
        end
        
        regr_coef = zeros([2 length(coordinates)]);
        parfor ii=1:length(coordinates)
            A = shifted_regr(ii,:);
            B = clean_voxel_ts(ii,:); B(isnan(A)) = []; A(isnan(A)) = [];
            C = [ones([length(A) 1]) A'];
            regr_coef(:,ii)= C\B';
        end
        
        cCVR(1,coordinates) = regr_coef(2,:); %extract slope
        cCVR(cCVR > 20) = 0; cCVR(cCVR < -20) = 0; %cleanup base CVR map
        cCVR = reshape(cCVR, [xx yy zz]);
        
        %calculate statistics
        
        SSE = zeros([1 size(wb_voxel_ts,1)]); SST = zeros([1 size(wb_voxel_ts,1)]); cT = zeros([1 size(wb_voxel_ts,1)]);
        parfor ii=1:size(wb_voxel_ts,1)
            A = shifted_regr(ii,:);
            X = clean_voxel_ts(ii,:); X(isnan(A)) = []; A(isnan(A)) = [];
            Y = regr_coef(2,ii)*A' + ones([1 length(A)])'.*regr_coef(1,ii);
            cT(1,ii) = regr_coef(2,ii)./nanstd(X-Y'); %(t = beta/STDEVr)
            SSE(1,ii) = (norm(X - Y'))^2;
            SST(1,ii) = (norm(X-mean(X)))^2;
        end
        
        R2 = 1 - SSE./SST;
        cR2 = zeros([1 xx*yy*zz]); cSSE = zeros([1 xx*yy*zz]);
        cR2(1, coordinates) = R2; cR2 = reshape(cR2, [xx yy zz]);
        cSSE(1, coordinates) = SSE; cSSE = reshape(cSSE, [xx yy zz]);
        cTstat(1, coordinates) = cT; cTstat = reshape(cTstat, [xx yy zz]);
        
        switch pp
            case 1 %entdidal
                %save lag-adjusted CVR map
                saveImageData(mask.*cCVR, opts.headers.map, opts.corrCVRdir,'cCVR_map.nii.gz', datatype);
                %save lag-adjusted stats data
                saveImageData(mask.*cR2, opts.headers.map, opts.corrCVRdir,'cR2_map.nii.gz', datatype);
                saveImageData(mask.*cSSE, opts.headers.map, opts.corrCVRdir,'cSSE_map.nii.gz', datatype);
                saveImageData(mask.*cTstat, opts.headers.map, opts.corrCVRdir,'cTstat_map.nii.gz', datatype);
                %save lagregressor maps
                lagregressor = zeros([xx*yy*zz dyn*opts.interp_factor]);
                lagregressor(coordinates,:) = shifted_regr;
                if opts.save_rts
                    saveImageData(lagregressor, opts.headers.ts, opts.corrCVRdir, 'lagregressor_map.nii.gz', datatype);
                end
                CVR(3,:,:,:) = mask.*cCVR;
                TSTAT(3,:,:,:) = mask.*cTstat;
                RC(3,:,:,:) = mask.*cR2;
                clear cCVR R2 cR2 cSSE lagregressor shifted_regr regr_coef SST SSE
            case 2 %effective probe
                %save lag-adjusted CVR map
                saveImageData(mask.*cCVR, opts.headers.map, opts.corrCVRdir,'cCVR_eff_map.nii.gz', datatype);
                %save lag-adjusted stats data
                saveImageData(mask.*cR2, opts.headers.map, opts.corrCVRdir,'cR2_eff_map.nii.gz', datatype);
                saveImageData(mask.*cSSE, opts.headers.map, opts.corrCVRdir,'cSSE_eff_map.nii.gz', datatype);
                saveImageData(mask.*cTstat, opts.headers.map, opts.corrCVRdir,'cTstat_eff_map.nii.gz', datatype);
                %save lagregressor maps
                lagregressor = zeros([xx*yy*zz dyn*opts.interp_factor]);
                lagregressor(coordinates,:) = shifted_regr;
                if opts.save_rts
                    saveImageData(lagregressor, opts.headers.ts, opts.corrCVRdir, 'lagregressor_eff_map.nii.gz', datatype);
                end
                CVR(4,:,:,:) = mask.*cCVR;
                TSTAT(4,:,:,:) = mask.*cTstat;
                RC(4,:,:,:) = mask.*cR2;
                clear cCVR R2 cR2 cSSE lagregressor shifted_regr regr_coef SST SSE
                
        end
        disp('Stimulus response data is saved')
    end
end
%% if CVR maps are generated using the input and effective probe, generate robust CVR map
if opts.cvr_maps && opts.eff_probe && opts.trace_corr
    disp('Calculating robust response maps')
    [MT,IT] = max(abs(TSTAT),[],1,'omitnan'); IT = squeeze(IT);
    [MR,IR] = max(abs(RC),[],1,'omitnan'); IR = squeeze(IR);
    CVR = reshape(CVR,4,xx*yy*zz);
    IT = reshape(IT,1,numel(IT));
    IR = reshape(IR,1,numel(IR));
    robustIT = zeros(size(IT));
    robustIR = zeros(size(IR));
    for ii=1:length(IT)
        tmp = IT(1,ii);
        robustIT(1,ii) = CVR(tmp,ii);
        tmp = IR(1,ii);
        robustIR(1,ii) = CVR(tmp,ii);
    end
    robustIT = reshape( robustIT,size(mask));
    robustIR = reshape( robustIR,size(mask));
    %save images
    if opts.robustTstat; saveImageData(mask.*robustIT, opts.headers.map, opts.corrCVRdir,'robustCVR_TSTAT.nii.gz', datatype); end
    if opts.robustR; saveImageData(mask.*robustIR, opts.headers.map, opts.corrCVRdir,'robustCVR_R.nii.gz', datatype); end
    clear TSTAT RC CVR robustIT robustIR
    disp('Calculating robust lag maps')

else
    clear TSTAT RC CVR
end


%% perform shifted regressor GLM

if opts.glm_model
     [~,lags] = xcorr(newprobe,newprobe,opts.highlag,'normalized');  %%%%%% !!!!!!!!!!
      idx = find(lags<=opts.adjlowlag | lags>=opts.adjhighlag);
      %remove lag values lower that what is possible
      lags(idx)=[];
        
    q = cputime;
    if opts.trace_corr; qq = [1 2]; else; qq = 1; end %regress with probe y/n
    for pp = qq
        switch pp
            case 1
                %setup regression  using optimized regressor
                disp('performing GLM based lag analysis using OPTIMIZED regressor')
                regr = newprobe;
            case 2
                %setup regression using end-tidal CO2
                disp('performing GLM based lag analysis using INPUT regressor')
                regr = probe;
        end
        %removes linearly dependent components
        
        for ii=1:size(np,1); np(ii,:) = np(ii,:)-rescale(np(ii,:)); end
        [norm_np,idx]=licols(np);
        %setup regression matrix
        for ii = 1:size(lags,2)
            corr_regr = circshift(regr,lags(ii));
            if lags(ii)> 1
                corr_regr(1:lags(ii),1) = regr(1,1);
            elseif lags(ii)< 1
                corr_regr(end-lags(ii):end,1) = regr(1,end);
            else
                %do nothing because index is zero = zero shift
            end
            regr_matrix(ii,:) = corr_regr;
        end
        %perform regression at all lag times
        regr_coef = zeros([size(lags,2) (size(norm_np,2)+2) length(coordinates)]);
        
        %run GLM
        parfor ii=1:size(regr_matrix,1)
            A = regr_matrix(ii,:);
            C = [ones([length(A(1,~isnan(A))) 1]) norm_np(~isnan(A),:) A(1,~isnan(A))']
            regr_coef(ii,:,:)= C\wb_voxel_ts(:,~isnan(A))';
        end
        
        %extract beta and compute lags
        estimatrix = squeeze(regr_coef(:,end,:)); %betas for last regressor
        if opts.uni
            [maximatrix maxindex] = max(estimatrix); %mx beta value defines lag
        else
            [maximatrix maxindex] = max(abs(estimatrix));
        end
        lagmatrix = lags(maxindex);
        GLM_Estimate = zeros([xx*yy*zz,1]); GLM_Estimate(coordinates,:) = maximatrix;
        GLM_lags = zeros([xx*yy*zz,1]); GLM_lags(coordinates,:) = lagmatrix;
        GLM_lags = reshape(GLM_lags,[xx yy zz]);
        GLM_Estimate = reshape(GLM_Estimate,[xx yy zz]);
        
        estimaTS = zeros([xx*yy*zz,length(lags)]);
        estimaTS(coordinates,:) = estimatrix';
        estimaTS = reshape(estimaTS, [xx yy zz length(lags)]);
        
        tmp = opts.TR*(lagmatrix/opts.interp_factor);
        tmpLag = zeros([xx*yy*zz,1]);  tmpLag(coordinates,:) = tmp; clear tmp
        tmpLag = reshape(tmpLag,[xx yy zz]);
        
        %calculate statistics
        
        SSE = zeros([1 size(wb_voxel_ts,1)]); SST = zeros([1 size(wb_voxel_ts,1)]);  T = zeros([1 size(wb_voxel_ts,1)]);
        cR2 = zeros([1 xx*yy*zz]); cSSE = zeros([1 xx*yy*zz]); cTstat = zeros([1 xx*yy*zz]);
        
        parfor ii=1:size(wb_voxel_ts,1)
            A = regr_matrix(maxindex(ii),:)
            X = wb_voxel_ts(ii,:); X(isnan(A)) = []; A(isnan(A)) = [];
            Y = regr_coef(maxindex(ii),end,ii)*A' + regr_coef(maxindex(ii),1,ii);
            cT(1,ii) = regr_coef(maxindex(ii),end,ii)./nanstd(X-Y'); %(t = beta/STDEVr)
            SSE(1,ii) = (norm(X - Y'))^2;
            SST(1,ii) = (norm(X-mean(X)))^2;
        end
        
        R2 = 1 - SSE./SST;
        cR2 = zeros([1 xx*yy*zz]); cSSE = zeros([1 xx*yy*zz]);
        cR2(1, coordinates) = R2; cR2 = reshape(cR2, [xx yy zz]);
        cSSE(1, coordinates) = SSE; cSSE = reshape(cSSE, [xx yy zz]);
        cTstat(1, coordinates) = cT; cTstat = reshape(cTstat, [xx yy zz]);
        
        %% plot(contributions)
        %output = input * coefficients
        %input = output / coefficients
        %coef = C\X --> X = C*coef
        
        nuis = zeros([size(norm_np,2) length(coordinates)]);
        intcp = zeros([1 length(coordinates)]);
        beta = zeros([1 length(coordinates)]);
        
        parfor ii=1:length(coordinates)
            nuis(:,ii) = squeeze(regr_coef(maxindex(ii),2:end-1,ii));
            intcp(1,ii) = squeeze(regr_coef(maxindex(ii),1,ii));
            beta(1,ii) = squeeze(regr_coef(maxindex(ii),end,ii));
        end
        
        nuis_TS = zeros([size(norm_np') length(coordinates)]);
        
        parfor ii=1:size(norm_np,2)
            nuis_TS(ii,:,:) =  np(:,ii)*nuis(ii,:);
        end
        combi_TS = squeeze(sum(nuis_TS,1));
        disp('check this')
        %all of the regressors used (i.e. last regressor)
        regr_TS = nan(size(wb_voxel_ts'));
        parfor ii=1:length(coordinates)
            regr_TS(:,ii) =  squeeze(regr_matrix(maxindex(ii),:))*beta(1,ii);
        end
        
        %plot regressors
        figure;
        subplot(4,1,1); plot(nanmean(wb_voxel_ts,1)); title('Original Data')
        subplot(4,1,2); plot(nanmean(regr_TS,2)); title('Main EV Component')
        subplot(4,1,3); plot(nanmean(wb_voxel_ts,1)'-nanmean(combi_TS,2)); title('Original data minus Nuissance Signal')
        subplot(4,1,4); plot(nanmean(wb_voxel_ts,1)'-nanmean(combi_TS,2)-nanmean(regr_TS,2)); title('Residual Data + Error Term')
        
        %save maps - shift lags to set minimum lag to zero (might not be
        %appropriate with pathology)
        
        switch pp
            case 1
                saveImageData(mask.*GLM_Estimate, opts.headers.map, opts.glmlagdir, 'optiReg_ES.nii.gz', datatype);
                saveImageData(mask.*tmpLag, opts.headers.map, opts.glmlagdir, 'optiReg_lags.nii.gz', datatype);
                saveImageData(mask.*GLM_lags, opts.headers.map,opts.glmlagdir, 'uncor_optiReg_lags.nii.gz', datatype);
                saveImageData(estimaTS, opts.headers.ts, opts.glmlagdir, 'optiReg_estimatrix.nii.gz', datatype);
                saveImageData(mask.*cR2, opts.headers.map, opts.glmlagdir, 'optiReg_R2.nii.gz', datatype);
                saveImageData(mask.*cSSE, opts.headers.map, opts.glmlagdir, 'optiReg_SSE.nii.gz', datatype);
                saveImageData(mask.*cTstat, opts.headers.map, opts.glmlagdir, 'optiReg_Tstat.nii.gz', datatype);
                saveas(gcf,[opts.glmlagdir,'regression_optiReg.fig']);
            case 2
                %save maps using CO2 regressor
                saveImageData(mask.*GLM_Estimate, opts.headers.map, opts.glmlagdir, 'inputReg_ES.nii.gz', datatype);
                saveImageData(mask.*tmpLag, opts.headers.map, opts.glmlagdir, 'inputReg_lags.nii.gz', datatype);
                saveImageData(mask.*GLM_lags, opts.headers.map, opts.glmlagdir, 'uncor_inputReg_lags.nii.gz', datatype);
                saveImageData(estimaTS, opts.headers.ts, opts.glmlagdir, 'inputReg_estimatrix.nii.gz', datatype);
                saveImageData(mask.*cR2, opts.headers.map, opts.glmlagdir, 'inputReg_R2.nii.gz', datatype);
                saveImageData(mask.*cSSE, opts.headers.map, opts.glmlagdir, 'inputReg_SSE.nii.gz', datatype);
                saveImageData(mask.*cTstat, opts.headers.map, opts.glmlagdir, 'inputReg_Tstat.nii.gz', datatype);
                saveas(gcf,[opts.glmlagdir,'regression_inputReg.fig']);
        end
        
        %%%%% generate lag-corrected CVR maps %%%%
        
        
        disp('Generating lag-adjusted CVR maps based on GLM analysis')
        cCVR = zeros([1 xx*yy*zz]);
        shifted_regr = NaN([length(coordinates) dyn*opts.interp_factor]);
        
        regr = probe;
        for ii = 1:length(coordinates)
            corr_regr = circshift(regr',maxindex(ii));
            if maxindex > 0
                corr_regr(1,1:maxindex(ii)) = NaN;
            elseif index < 0
                corr_regr(1,end-maxindex:end) = NaN;
            else
                %do nothing because index is zero = zero shift
            end
            shifted_regr(ii,:) = corr_regr;
        end
        
        regr_coef = zeros([2 length(coordinates)]);
        parfor ii=1:length(coordinates)
            A = shifted_regr(ii,:);
            B = clean_voxel_ts(ii,:); B(isnan(A)) = []; A(isnan(A)) = [];
            C = [ones([length(A) 1]) A'];
            regr_coef(:,ii)= C\B';
        end
        
        cCVR(1,coordinates) = regr_coef(end,:); %extract slope
        cCVR(cCVR > 20) = 0; cCVR(cCVR < -20) = 0; %cleanup base CVR map
        cCVR = reshape(cCVR, [xx yy zz]);
        %save lag-adjusted CVR map
        if pp == 1
            saveImageData(mask.*cCVR, opts.headers.map, opts.glmCVRdir, 'optiReg_cCVR.nii.gz', datatype);
        else
            saveImageData(mask.*cCVR, opts.headers.map, opts.glmCVRdir, 'inputReg_cCVR.nii.gz', datatype);
        end
        
        %calculate statistics
        SSE = zeros([1 size(wb_voxel_ts,1)]); SST = zeros([1 size(wb_voxel_ts,1)]); cT = zeros([1 size(wb_voxel_ts,1)]);
        parfor ii=1:size(wb_voxel_ts,1)
            A = shifted_regr(ii,:);
            X = clean_voxel_ts(ii,:); X(isnan(A)) = []; A(isnan(A)) = [];
            Y = regr_coef(2,ii)*A' + regr_coef(1,ii);
            cT(1,ii) = regr_coef(2,ii)./nanstd(X-Y'); %(t = beta/STDEVr)
            SSE(1,ii) = (norm(X - Y'))^2;
            SST(1,ii) = (norm(X-mean(X)))^2;
        end
        
        R2 = 1 - SSE./SST;
        cR2 = zeros([1 xx*yy*zz]); cSSE = zeros([1 xx*yy*zz]);  cTstat = zeros([1 xx*yy*zz]);
        cR2(1, coordinates) = R2'; cR2 = reshape(cR2, [xx yy zz]);
        cSSE(1, coordinates) = SSE; cSSE = reshape(cSSE, [xx yy zz]);
        cTstat(1, coordinates) = cT; cTstat = reshape(cTstat, [xx yy zz]);
        %save lag-adjusted stats data
        if pp == 1
            saveImageData(mask.*cR2, opts.headers.map, opts.glmCVRdir, 'optiReg_cR2.nii.gz', datatype);
            saveImageData(mask.*cSSE, opts.headers.map, opts.glmCVRdir, 'optiReg_cSSE.nii.gz', datatype);
            saveImageData(mask.*cTstat, opts.headers.map, opts.glmCVRdir,'optiReg_cTstat.nii.gz', datatype);
        else
            saveImageData(mask.*cR2, opts.headers.map, opts.glmCVRdir, 'inputReg_cR2.nii.gz', datatype);
            saveImageData(mask.*cSSE, opts.headers.map, opts.glmCVRdir, 'inputReg_cSSE.nii.gz', datatype);
            saveImageData(mask.*cTstat, opts.headers.map, opts.glmCVRdir,'inputReg_cTstat.nii.gz', datatype);
        end
        
    end
    disp(['finished running GLM analysis in: ',int2str((cputime-q)/60),' minutes'])
    
end

