% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% sections of this code were contributed by Allen A. Champagne, a.champagne@queensu.ca
% <lagCVR: calculates hemodynamic parameter maps with associated statistical maps >
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

function [newprobe, maps] = lagCVR(refmask, mask, data, probe, nuisance, opts)
% This function calculates hemodynamic parameter maps with associated
% statistical maps. For full usage details download the user manual from
% https://www.seevr.nl/download-seevr/ or look at usage tutorials
% https://www.seevr.nl/tutorials/
refmask = logical(refmask); mask = logical(mask);
probe = double(probe);
data = double(data);

warning('off');
global opts;
data(isnan(data)) = 0; data(isinf(data)) = 0;

%check for stats toolbox
if license('test','Statistics_toolbox') == 0; opts.pca_analysis = 0; end

maps = struct();
if iscolumn(probe); else; probe = probe'; end

if isempty(nuisance); np = []; else
    test1 = nuisance(1,:); test2 = nuisance(:,1);
    if length(test1) > length(test2); nuisance = nuisance'; end
    for ii=1:size(nuisance,2); np(:,ii) = rescale(interp(nuisance(:,ii),opts.interp_factor),-1,1); end
    %removes linearly dependent components
    %norm_np = np;
    [norm_np,idx]=licols(np);
    clear test1 test2
end

%setup default parameters
if isfield(opts,'gpu'); else; opts.gpu = 0; end                           %use gpu
if isfield(opts,'niiwrite'); else; opts.niiwrite = 0; end                 %depending on how data is loaded this can be set to 1 to use native load/save functions
if isfield(opts,'plot'); else; opts.plot = 0; end                         %show or hide selected plots
if isfield(opts,'prewhite'); else; opts.prewhite = 0; end                 %zero mean and unit variance of data. normally leads to bad results
if isfield(opts,'interp_factor'); else; opts.interp_factor = 4; end       %factor by which to temporally interpolate data. Better for picking up lags between TR
if isfield(opts,'load_probe'); else; opts.load_probe = 0; end             %saves time for creating regressor if one with the correct length already exists
if isfield(opts,'save_rts'); else; opts.save_rts = 0; end                 %save correlation timeseries - can be used to visualize lags
if isfield(opts,'rescale_probe'); else; opts.rescale_probe = 1; end       %rescaling may be helful for initial refinement run
if isfield(opts,'trace_corr'); else; opts.trace_corr = 1; end             %perform additional correlation with gas trace (on top of optimized regressor)
if isfield(opts,'refine_regressor'); else; opts.refine_regressor = 1; end %refine BOLD regressor. If '0' a straight correlation will be done with probe
if isfield(opts,'pca_analysis'); else; opts.pca_analysis = 1; end         %PCA analysis to optimize BOLD regressor
if isfield(opts,'corr_model'); else; opts.corr_model = 1; end             %perform correlation based analysis
if isfield(opts,'cvr_maps'); else; opts.cvr_maps = 1; end                 %generate CVR maps based on regression of BOLD signal with gas probe
if isfield(opts,'eff_probe'); else; opts.eff_probe = 0; end               %uses effective probe to calculate probe response
if isfield(opts,'glm_model'); else; opts.glm_model = 0; end               %perform GLM based lag regression. This works well for high temporal resolution data
if isfield(opts,'uni'); else; opts.uni = 0; end                           %if this is set to 1 negative correlations will be ignored
if isfield(opts,'norm_regr'); else; opts.norm_regr = 0; end               %normalize the regressor for correlation analysis between 0-1
if isfield(opts,'robust'); else; opts.robust = 0; end                     %calculate robust lag and CVR
if opts.niiwrite
    if isfield(opts.info,'rts'); else; opts.info.rts = opts.info.ts; end
end

%important parameters for results
if isfield(opts,'corr_thresh'); else; opts.corr_thresh = 0.7; end         %threshold by which to accept correlated voxels during the refinement of the BOLD regressor
%refinement
if isfield(opts,'lowerlagthresh'); else; opts.lowerlagthresh = -2; end    %lower threshold for refinement (generaly -2 to 0)
if isfield(opts,'upperlagthresh'); else; opts.upperlagthresh = 2; end     %upper threshold for refinement (generaly 0 to +2)

%account for interpolation factor
opts.adjlowerthresh = opts.lowerlagthresh*opts.interp_factor;
opts.adjupperthresh = opts.upperlagthresh*opts.interp_factor;

%correlation
if isfield(opts,'lowlag'); else; opts.lowlag = -3; end                    %lower threshold for correlation (generaly -3 to 0)
if isfield(opts,'highlag'); else; opts.highlag = 60; end                  %upper threshold for correlation (in healthy up to 20-40, in disease 60-90)

%account for interpolation factor
opts.adjlowlag = opts.lowlag*opts.interp_factor; %setup lower lag limit; negative for misalignment and noisy correlation
opts.adjhighlag = opts.highlag*opts.interp_factor; %setups upper lag limit; allow for long lags associated with pathology

%check
if opts.load_probe == 0 && opts.refine_regressor == 0 && opts.trace_corr == 0
    disp('check options; stopping analysis')
    opts.glm_model = 0;
    opts.corr_model = 0;
    opts.cvr_maps = 0;
end
if opts.trace_corr == 0 && opts.robust == 1
    disp('For robust analysis, set opts.trace_corr = 1')
    disp('...continuing without creating robust maps')
    opts.robust = 0;
end

%setup save directories
if isfield(opts,'resultsdir'); else; opts.resultsdir = fullfile(pwd); end
if opts.corr_model; opts.corrlagdir = fullfile(opts.resultsdir,'corrLAG'); mkdir(opts.corrlagdir); end
if opts.corr_model && opts.cvr_maps; opts.corrCVRdir = fullfile(opts.corrlagdir,'CVR'); mkdir(opts.corrCVRdir); end
if opts.corr_model == 0; opts.cvr_maps = 0; end
if opts.glm_model
    opts.glmlagdir = fullfile(opts.resultsdir,'glmLAG'); mkdir(opts.glmlagdir);
    if opts.cvr_maps; opts.glmCVRdir = fullfile(opts.glmlagdir,'CVR'); mkdir(opts.glmCVRdir); end
end

%optional image outputs
if isfield(opts, 'robustTstat'); else; opts.robustTstat = 1; end
if isfield(opts, 'robustR'); else; opts.robustR = 0; end

cd(opts.resultsdir);
datatype = 16;
[xx, yy, zz, dyn] = size(data);
orig_regr = probe;

t = cputime;

%% grab coordinates
% WB coordinates
[orig_voxel_ts, coordinates] = grabTimeseries(data, mask);
% GM coordinates
[gm_voxel_ts, gmcoordinates] = grabTimeseries(data, refmask);

if opts.prewhite
    gm_voxel_ts = gm_voxel_ts'; parfor ii=1:size(gmcoordinates,1); [gm_voxel_ts(:,ii), ~, ~] = prewhiten(gm_voxel_ts(:,ii)); end; gm_voxel_ts = gm_voxel_ts';
    pw_voxel_ts = orig_voxel_ts'; parfor ii=1:size(coordinates,1); [pw_voxel_ts(:,ii), ~, ~] = prewhiten(pw_voxel_ts(:,ii)); end; pw_voxel_ts = pw_voxel_ts';
    [probe, ~, ~] = prewhiten(probe);
end

%interpolate variables
probe = interp(probe,opts.interp_factor); %interpolate data by a factor of the number of slices for more accurate timing
orig_regr = interp(orig_regr,opts.interp_factor); %save input regressor to generate corrected CVR maps
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
        
        %%
        
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
            
            a=gm_voxel_ts_nonan;
            a2=mat2cell(a,ones(1,size(a,1)),size(a,2)); %turn into cells so treat all rows independently
            
            % the order here matters because of the SIGN of the lag
            % this way we slide the probe to the ts
            b2=cellfun(@(a2) xcorr(a2,newprobe,'coeff'),a2,'uniformoutput',false);
            corr_probe=cell2mat(b2); %reformat into correlation matrix
            corr_probe = corr_probe' ;
            
            %remove low and high lags ignoring unreasonable and long lag times
            idx = find(lags<=opts.adjlowerthresh | lags>=opts.adjupperthresh);
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
            checklag(ts_lag>=opts.adjlowerthresh & ts_lag<=opts.adjupperthresh)=0;
            checklag(ts_lag<opts.adjlowerthresh | ts_lag>opts.adjupperthresh)=1;
            
            resultant_filter = zeros(size(checklag));
            resultant_filter(checklag==1 | checkcorr==1)=1;
            remove = find(resultant_filter==1);
            
            % removing timeseries and lags outside of range
            gm_voxel_ts_nonan_filt=gm_voxel_ts_nonan';
            gm_voxel_ts_nonan_filt(:,remove)=[];
            ts_lag_filt=ts_lag;
            ts_lag_filt(remove)=[];
            
            % re-aligning all timeseries based on lag for PCA
            tic
            new_shifted_ts = zeros(length(ts_lag_filt), length(probe));
            for hh = 1:length(ts_lag_filt)
                new_shifted_ts(hh,:) = circshift(gm_voxel_ts_nonan_filt(:,hh),ts_lag_filt(hh));
                if ts_lag_filt(hh) > 0
                    new_shifted_ts(hh,1:ts_lag_filt(hh)) = gm_voxel_ts_nonan_filt(1,1);
                elseif ts_lag_filt(hh) < 0
                    new_shifted_ts(hh,end-ts_lag_filt(hh):end) = gm_voxel_ts_nonan_filt(1,end);
                else
                    new_shifted_ts(hh,:) = gm_voxel_ts_nonan_filt(:,hh);
                end
            end
            toc
            
            disp('finished correlation, estimating new probe')
            
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
                    disp(['the opts.corr_thresh parameter is currently set at: ',num2str(opts.corr_thresh)])
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
        for ii=1:size(keep_probes,1); legendInfo{ii} = ['probe ',int2str(ii)]; end
        subplot(3,1,1);
        plot(probe,'LineWidth',2); title('input probe' )
        subplot(3,1,2);
        plot(keep_probes','LineWidth',2); title('probe iterations');
        legend(legendInfo)
        subplot(3,1,3);
        plot(newprobe,'LineWidth',2); title('optimized probe');
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
    idx = find(lags<=opts.adjlowlag | lags>=opts.adjlowlag);
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
        
        b2=cellfun(@(a2) xcorr(a2,regr,opts.adjhighlag,'coeff'),a2,'uniformoutput',false);
        corr_probe=cell2mat(b2); %reformat into correlation matrix
        corr_probe = corr_probe';
        %use highlag to restrict the number of correlations that need to be done
        [~,lags] = xcorr(gm_voxel_ts_nonan(1,:),regr,opts.adjhighlag,'coeff');  %%%%%% !!!!!!!!!!
        
        %save correlations over time
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
        tmpLag = opts.TR*(lag_map/opts.interp_factor); tmpLag(1,coordinates) = tmpLag(1,coordinates) + abs(min(tmpLag(:))); %puts lag maps back into seconds
        lag_map = reshape(lag_map,[xx yy zz]);
        tmpLag = reshape(tmpLag,[xx yy zz]);
        %save lag and r maps
        
        switch pp
            case 1
                if opts.niiwrite
                    cd(opts.corrlagdir);
                    niftiwrite(cast(mask.*tmpLag,opts.mapDatatype),'lag_map',opts.info.map);
                    niftiwrite(cast(mask.*lag_map,opts.mapDatatype),'uncorr_lag_map',opts.info.map);
                    niftiwrite(cast(mask.*r_map,opts.mapDatatype),'r_map',opts.info.map);
                else
                    saveImageData(mask.*tmpLag, opts.headers.map, opts.corrlagdir, 'lag_map.nii.gz', datatype);
                    saveImageData(mask.*lag_map, opts.headers.map, opts.corrlagdir, 'uncorr_lag_map.nii.gz', datatype);
                    saveImageData(mask.*r_map, opts.headers.map, opts.corrlagdir, 'r_map.nii.gz', datatype);
                end
                maps.XCORR.lag_opti = tmpLag;
                maps.XCORR.uncorrlag_opti = lag_map;
                maps.XCORR.r_opti = r_map;
                
                if opts.save_rts
                    if opts.niiwrite
                        cd(opts.corrlagdir);
                        disp('saving correlation timeseries based on optimized regressor')
                        niftiwrite(cast(corr_ts,opts.mapDatatype),'r_ts',opts.info.rts);
                        clear Trcorr_ts
                    else
                        disp('saving correlation timeseries based on optimized regressor')
                        saveImageData(corr_ts, headers.ts, opts.corrlagdir,  'r_ts.nii.gz', datatype)
                        clear Trcorr_ts
                    end
                end
                if opts.robust
                    LAG(1,:,:,:) = mask.*mask.*tmpLag;
                    RL(1,:,:,:) = mask.*r_map;
                end
                clear tmpLag lag_map r_map
            case 2
                if opts.niiwrite
                    cd(opts.corrlagdir);
                    niftiwrite(cast(mask.*tmpLag,opts.mapDatatype),'lag_map_probe',opts.info.map);
                    niftiwrite(cast(mask.*lag_map,opts.mapDatatype),'uncorr_lag_map_probe',opts.info.map);
                    niftiwrite(cast(mask.*r_map,opts.mapDatatype),'r_map_probe',opts.info.map);
                else
                    saveImageData(mask.*tmpLag, opts.headers.map, opts.corrlagdir,  'lag_map_probe.nii.gz', datatype);
                    saveImageData(mask.*lag_map, opts.headers.map,opts.corrlagdir, 'uncorr_lag_map_probe.nii.gz', datatype);
                    saveImageData(mask.*r_map, opts.headers.map, opts.corrlagdir, 'r_map_probe.nii.gz', datatype);
                end
                maps.XCORR.lag_input = tmpLag;
                maps.XCORR.uncorrlag_input = lag_map;
                maps.XCORR.r_input = r_map;
                
                if opts.save_rts
                    if opts.niiwrite
                        cd(opts.corrlagdir);
                        disp('saving correlation timeseries based on optimized regressor')
                        niftiwrite(cast(corr_ts,opts.mapDatatype),'r_ts_probe',opts.info.rts);
                        clear Trcorr_ts
                    else
                        disp('saving correlation timeseries based on input probe')
                        saveImageData(Trcorr_ts, opts.headers.ts, opts.corrlagdir, 'r_ts_probe.nii.gz', datatype);
                        clear Trcorr_ts
                    end
                end
                if opts.robust
                    LAG(2,:,:,:) = mask.*mask.*tmpLag;
                    RL(2,:,:,:) = mask.*r_map;
                end
                clear tmpLag lag_map r_map
        end
        clear Trcorr_ts orig_voxel_ts
    end
    
    disp(['Lag, regressor and r_maps were created in: ',int2str((cputime-t)/60),' minutes'])
    
    clear gm_voxel_ts_nonan m corr_ts regr
end


%% if correlating both input regressor and optimized probe
if opts.trace_corr && opts.robust
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
    if opts.niiwrite
        cd(opts.corrlagdir);
        niftiwrite(cast(mask,opts.mapDatatype).*robustIR,'robustLAG_r',opts.info.map);
    else
        saveImageData(mask.*robustIR, opts.headers.map, opts.corrlagdir,'robustLAG_r.nii.gz', datatype);
    end
    maps.XCORR.robust_lag_r = robustIR;
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
                regr = orig_regr;
            case 2
                %setup regression using end-tidal CO2
                rs_newprobe = rescale(newprobe, 0.0001,1 );
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
        regr_coef = C\clean_voxel_ts';
        
        bCVR(1,coordinates) = regr_coef(2,:); %extract slope
        bCVR(bCVR > 10) = 0; bCVR(bCVR < -10) = 0; %cleanup base CVR map
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
                if opts.niiwrite
                    cd(opts.corrCVRdir);
                    niftiwrite(cast(mask.*bCVR,opts.mapDatatype),'bCVR_map',opts.info.map);
                    niftiwrite(cast(mask.*bR2,opts.mapDatatype),'bR2_map',opts.info.map);
                    niftiwrite(cast(mask.*bSSE,opts.mapDatatype),'bSSE2_map',opts.info.map);
                    niftiwrite(cast(mask.*bTstat,opts.mapDatatype),'bTstat_map',opts.info.map);
                else
                    saveImageData(mask.*bCVR, opts.headers.map, opts.corrCVRdir, 'bCVR_map.nii.gz', datatype);
                    saveImageData(mask.*bR2, opts.headers.map, opts.corrCVRdir,'bR2_map.nii.gz', datatype);
                    saveImageData(mask.*bSSE, opts.headers.map, opts.corrCVRdir,'bSSE_map.nii.gz', datatype);
                    saveImageData(mask.*bTstat, opts.headers.map, opts.corrCVRdir,'bTstat_map.nii.gz', datatype);
                end
                maps.XCORR.CVR.bCVR = bCVR;
                maps.XCORR.CVR.bR2 = bR2;
                maps.XCORR.CVR.bSSE = bSSE;
                maps.XCORR.CVR.bTstat = bSSE;
                
                disp('Base CVR, r2, Tstat and SSE maps were created using entidal regressor')
                CVR(1,:,:,:) = mask.*bCVR;
                TSTAT(1,:,:,:) = mask.*bTstat;
                RC(1,:,:,:) = mask.*bR2;
                clear bCVR bR2 bT bTstat bSSE SSE SST R2 X Y r regr_coef A C s
                
            case 2 %effective probe
                if opts.niiwrite
                    opts.corrCVRdir
                    niftiwrite(cast(mask.*bCVR,opts.mapDatatype),'bCVR_eff_map',opts.info.map);
                    niftiwrite(cast(mask.*bR2,opts.mapDatatype),'bR2_eff_map',opts.info.map);
                    niftiwrite(cast(mask.*bSSE,opts.mapDatatype),'bSSE2_eff_map',opts.info.map);
                    niftiwrite(cast(mask.*bTstat,opts.mapDatatype),'bTstat_eff_map',opts.info.map);
                else
                    saveImageData(mask.*bCVR, opts.headers.map, opts.corrCVRdir, 'bCVR_eff_map.nii.gz', datatype);
                    saveImageData(mask.*bR2, opts.headers.map, opts.corrCVRdir,'bR2_eff_map.nii.gz', datatype);
                    saveImageData(mask.*bSSE, opts.headers.map, opts.corrCVRdir,'bSSE_eff_map.nii.gz', datatype);
                    saveImageData(mask.*bTstat, opts.headers.map, opts.corrCVRdir,'bTstat_eff_map.nii.gz', datatype);
                end
                maps.XCORR.CVR.bCVR_eff = bCVR;
                maps.XCORR.CVR.bR2_eff = bR2;
                maps.XCORR.CVR.bSSE_eff = bSSE;
                maps.XCORR.CVR.bTstat_eff = bSSE;
                
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
                corr_regr(1,1:index(ii)) = regr(1,1); %formerly NaN
            elseif index < 0
                corr_regr(1,end-index:end) = regr(1,end); %formerly NaN
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
        cCVR(cCVR > 10) = 0; cCVR(cCVR < -10) = 0; %cleanup base CVR map
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
                if opts.niiwrite
                    cd(opts.corrCVRdir);
                    niftiwrite(cast(mask,opts.mapDatatype).*cCVR,'cCVR_map',opts.info.map);
                    niftiwrite(cast(mask,opts.mapDatatype).*cR2,'cR2_map',opts.info.map);
                    niftiwrite(cast(mask,opts.mapDatatype).*cSSE,'cSSE2_map',opts.info.map);
                    niftiwrite(cast(mask,opts.mapDatatype).*cTstat,'cTstat_map',opts.info.map);
                else
                    saveImageData(mask.*cCVR, opts.headers.map, opts.corrCVRdir,'cCVR_map.nii.gz', datatype);
                    saveImageData(mask.*cR2, opts.headers.map, opts.corrCVRdir,'cR2_map.nii.gz', datatype);
                    saveImageData(mask.*cSSE, opts.headers.map, opts.corrCVRdir,'cSSE_map.nii.gz', datatype);
                    saveImageData(mask.*cTstat, opts.headers.map, opts.corrCVRdir,'cTstat_map.nii.gz', datatype);
                end
                maps.XCORR.CVR.cCVR = cCVR;
                maps.XCORR.CVR.cR2 = cR2;
                maps.XCORR.CVR.cSSE = cSSE;
                maps.XCORR.CVR.cTstat = cTstat;
                
                %save lagregressor maps
                lagregressor = zeros([xx*yy*zz dyn*opts.interp_factor]);
                lagregressor(coordinates,:) = shifted_regr;
                if opts.save_rts
                    if opts.niiwrite
                        cd(opts.corrCVRdir);
                        niftiwrite(cast(lagregressor,opts.mapDatatype),'lagregressor_map',opts.info.rts);
                    else
                        saveImageData(lagregressor, opts.headers.ts, opts.corrCVRdir, 'lagregressor_map.nii.gz', datatype);
                    end
                end
                if opts.robust
                    CVR(3,:,:,:) = mask.*cCVR;
                    TSTAT(3,:,:,:) = mask.*cTstat;
                    RC(3,:,:,:) = mask.*cR2;
                end
                clear cCVR R2 cR2 cSSE lagregressor shifted_regr regr_coef SST SSE
            case 2 %effective probe
                if opts.niiwrite
                    cd(opts.corrCVRdir);
                    niftiwrite(cast(mask.*cCVR,opts.mapDatatype),'cCVR_eff_map',opts.info.map);
                    niftiwrite(cast(mask.*cR2,opts.mapDatatype),'cR2_eff_map',opts.info.map);
                    niftiwrite(cast(mask.*cSSE,opts.mapDatatype),'cSSE2_eff_map',opts.info.map);
                    niftiwrite(cast(mask.*cTstat,opts.mapDatatype),'cTstat_eff_map',opts.info.map);
                else
                    saveImageData(mask.*cCVR, opts.headers.map, opts.corrCVRdir,'cCVR_eff_map.nii.gz', datatype);
                    saveImageData(mask.*cR2, opts.headers.map, opts.corrCVRdir,'cR2_eff_map.nii.gz', datatype);
                    saveImageData(mask.*cSSE, opts.headers.map, opts.corrCVRdir,'cSSE_eff_map.nii.gz', datatype);
                    saveImageData(mask.*cTstat, opts.headers.map, opts.corrCVRdir,'cTstat_eff_map.nii.gz', datatype);
                end
                maps.XCORR.CVR.cCVR_eff = cCVR;
                maps.XCORR.CVR.cR2_eff = cR2;
                maps.XCORR.CVR.cSSE_eff = cSSE;
                maps.XCORR.CVR.cTstat_eff = cTstat;
                
                
                %save lagregressor maps
                lagregressor = zeros([xx*yy*zz dyn*opts.interp_factor]);
                lagregressor(coordinates,:) = shifted_regr;
                if opts.save_rts
                    if opts.niiwrite
                        cd(opts.corrCVRdir);
                        niftiwrite(cast(lagregressor,opts.mapDatatype),'lagregressor_map_eff',opts.info.rts);
                    else
                        saveImageData(lagregressor, opts.headers.ts, opts.corrCVRdir, 'lagregressor_eff_map.nii.gz', datatype);
                    end
                end
                if opts.robust
                    CVR(4,:,:,:) = mask.*cCVR;
                    TSTAT(4,:,:,:) = mask.*cTstat;
                    RC(4,:,:,:) = mask.*cR2;
                end
                clear cCVR R2 cR2 cSSE lagregressor shifted_regr regr_coef SST SSE
        end
        disp('Stimulus response data is saved')
    end
end
%% if CVR maps are generated using the input and effective probe, generate robust CVR map
if opts.cvr_maps && opts.eff_probe && opts.trace_corr && opts.robust
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
    if opts.robustTstat
        if opts.niiwrite
            cd(opts.corrCVRdir);
            niftiwrite(cast(mask.*robustIT,opts.mapDatatype),'robustCVR_TSTAT',opts.info.map);
        else
            saveImageData(mask.*robustIT, opts.headers.map, opts.corrCVRdir,'robustCVR_TSTAT.nii.gz', datatype);
        end
        maps.XCORR.CVR.robustCVR_TSTAT = robustIT;
    end
    if opts.robustR
        if opts.niiwrite
            cd(opts.corrCVRdir);
            niftiwrite(cast(mask.*robustIR,opts.mapDatatype),'robustCVR_R',opts.info.map);
        else
            saveImageData(mask.*robustIR, opts.headers.map, opts.corrCVRdir,'robustCVR_R.nii.gz', datatype);
        end
        maps.XCORR.CVR.robustCVR_R = robustIR;
    end
    clear TSTAT RC CVR robustIT robustIR
    disp('Calculating robust lag maps')
    
else
    clear TSTAT RC CVR
end


%% perform shifted regressor GLM

if opts.glm_model
    [~,lags] = xcorr(newprobe,newprobe,opts.adjhighlag,'coeff');  %%%%%% !!!!!!!!!!
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
        
        %setup regression matrix
        for ii = 1:size(lags,2)
            corr_regr = circshift(regr,lags(ii));
            if lags(ii)> 1
                corr_regr(1:lags(ii),1) = regr(1,1); %formerly NaN
            elseif lags(ii)< 1
                corr_regr(end-lags(ii):end,1) = regr(1,end); %formerly NaN
            else
                %do nothing because index is zero = zero shift
            end
            regr_matrix(ii,:) = corr_regr;
        end
        %perform regression at all lag times
        
        if isempty(np) || nnz(np) == 0
            if opts.gpu
                regr_coef = zeros([size(lags,2) 2 length(coordinates)]);
                wb_voxel_ts = gpuArray(wb_voxel_ts);
                for ii=1:size(regr_matrix,1)
                    A = gpuArray(regr_matrix(ii,:));
                    C = gpuArray([ones([length(A(1,~isnan(A))) 1]) A(1,~isnan(A))']);
                    regr_coef(ii,:,:)= gather(C\wb_voxel_ts(:,~isnan(A))');
                end
                wb_voxel_ts = gather(wb_voxel_ts);
            else
                regr_coef = zeros([size(lags,2) 2 length(coordinates)]);
                parfor ii=1:size(regr_matrix,1)
                    A = regr_matrix(ii,:);
                    C = [ones([length(A(1,~isnan(A))) 1]) A(1,~isnan(A))'];
                    regr_coef(ii,:,:)= C\wb_voxel_ts(:,~isnan(A))';
                end
            end
            %reconstruct signals to calculate R2
            
            R2 = zeros([1 size(regr_matrix,1)]);
            maxindex = zeros([1 size(wb_voxel_ts,1)]);
            beta = zeros([1 size(wb_voxel_ts,1)]);
            rsquared = zeros([1 size(wb_voxel_ts,1)]);
            
            for ii = 1:size(wb_voxel_ts,1)
                tcoef = regr_coef(:,:,ii);
                X =  wb_voxel_ts(ii,:);
                for jj=1:size(regr_matrix,1)
                    tY = tcoef(jj,2).*regr_matrix(jj,:) + tcoef(jj,1);
                    R2(jj) = corr2(tY,X)^2;
                end
                [M, I] = max(R2');
                maxindex(1,ii) = I;
                rsquared(1,ii) = M;
                beta(1,ii) = regr_coef(I,2,ii);
            end
            beta(beta > 10) = 0; beta(beta < -10) = 0;
        else
            regr_coef = zeros([size(lags,2) (size(norm_np,2)+2) length(coordinates)]);
            %run GLM
            if opts.gpu
                wb_voxel_ts = gpuArray(wb_voxel_ts);
                for ii=1:size(regr_matrix,1)
                    A = gpuArray(regr_matrix(ii,:));
                    C = gpuArray([ones([length(A(1,~isnan(A))) 1]) norm_np(~isnan(A),:) A(1,~isnan(A))']);
                    regr_coef(ii,:,:)= gather(C\wb_voxel_ts(:,~isnan(A))');
                end
                wb_voxel_ts = gather(wb_voxel_ts);
            else
                parfor ii=1:size(regr_matrix,1)
                    A = regr_matrix(ii,:);
                    C = [ones([length(A(1,~isnan(A))) 1]) norm_np(~isnan(A),:) A(1,~isnan(A))'];
                    regr_coef(ii,:,:)= C\wb_voxel_ts(:,~isnan(A))';
                end
            end
            %reconstruct signals to calculate R2
            
            R2 = zeros([1 size(regr_matrix,1)]);
            maxindex = zeros([1 size(wb_voxel_ts,1)]);
            beta = zeros([1 size(wb_voxel_ts,1)]);
            rsquared = zeros([1 size(wb_voxel_ts,1)]);
            nuis = zeros([size(norm_np,1) 1]);
            
            for ii = 1:size(wb_voxel_ts,1)
                tcoef = regr_coef(:,:,ii);
                X =  wb_voxel_ts(ii,:);
                for jj=1:size(regr_matrix,1)
                    tnuis = norm_np*squeeze(tcoef(jj,2:end-1)'); % check
                    tY =  tcoef(jj,end).*regr_matrix(jj,:)+ (nuis + sum(tnuis,2))' + tcoef(jj,1);
                    R2(jj) = corr2(tY,X)^2;
                end
                [M, I] = max(R2');
                maxindex(1,ii) = I;
                rsquared(1,ii) = M;
                beta(1,ii) = regr_coef(I,2,ii);
            end
            beta(beta > 30) = 0; beta(beta < -30) = 0;
        end
        
        lagmatrix = lags(maxindex);
        GLM_Estimate = zeros([xx*yy*zz,1]); GLM_Estimate(coordinates,:) = beta;
        GLM_lags = zeros([xx*yy*zz,1]); GLM_lags(coordinates,:) = lagmatrix;
        GLM_lags = reshape(GLM_lags,[xx yy zz]);
        GLM_Estimate = reshape(GLM_Estimate,[xx yy zz]);
        
        tmp = opts.TR*(lagmatrix/opts.interp_factor);
        tmpLag = zeros([xx*yy*zz,1]);  tmpLag(coordinates,:) = tmp + abs(min(tmp(:))); clear tmp
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
        
        
        switch pp
            case 1
                if opts.niiwrite
                    cd(opts.glmlagdir);
                    niftiwrite(cast(mask.*GLM_Estimate,opts.mapDatatype),'optiReg_beta',opts.info.map);
                    niftiwrite(cast(mask.*tmpLag,opts.mapDatatype),'optiReg_lags',opts.info.map);
                    niftiwrite(cast(mask.*GLM_lags,opts.mapDatatype),'uncor_optiReg_lags',opts.info.map);
                    niftiwrite(cast(mask.*cR2,opts.mapDatatype),'optiReg_R2',opts.info.map);
                    niftiwrite(cast(mask.*cSSE,opts.mapDatatype),'optiReg_SSE',opts.info.map);
                    niftiwrite(cast(mask.*cTstat,opts.mapDatatype),'optiReg_Tstat',opts.info.map);
                else
                    saveImageData(mask.*GLM_Estimate, opts.headers.map, opts.glmlagdir, 'optiReg_beta.nii.gz', datatype);
                    saveImageData(mask.*tmpLag, opts.headers.map, opts.glmlagdir, 'optiReg_lags.nii.gz', datatype);
                    saveImageData(mask.*GLM_lags, opts.headers.map,opts.glmlagdir, 'uncor_optiReg_lags.nii.gz', datatype);
                    saveImageData(mask.*cR2, opts.headers.map, opts.glmlagdir, 'optiReg_R2.nii.gz', datatype);
                    saveImageData(mask.*cSSE, opts.headers.map, opts.glmlagdir, 'optiReg_SSE.nii.gz', datatype);
                    saveImageData(mask.*cTstat, opts.headers.map, opts.glmlagdir, 'optiReg_Tstat.nii.gz', datatype);
                end
                maps.GLM.optiReg_ES = GLM_Estimate;
                maps.GLM.optiReg_lags = tmpLag;
                maps.GLM.uncor_optiReg_lags = GLM_lags;
                maps.GLM.optiReg_R2 = cR2;
                maps.GLM.optiReg_SSE = cSSE;
                maps.GLM.optiReg_Tstat = cTstat;
                %if opts.plot; saveas(gcf,[opts.glmlagdir,'regression_optiReg.fig']); end
                
            case 2 %save maps using CO2 regressor
                if opts.niiwrite
                    cd(opts.glmlagdir);
                    niftiwrite(cast(mask.*GLM_Estimate,opts.mapDatatype),'inputReg_beta',opts.info.map);
                    niftiwrite(cast(mask.*tmpLag,opts.mapDatatype),'inputReg_lags',opts.info.map);
                    niftiwrite(cast(mask.*GLM_lags,opts.mapDatatype),'uncor_inputReg_lags',opts.info.map);
                    niftiwrite(cast(mask.*cR2,opts.mapDatatype),'inputReg_R2',opts.info.map);
                    niftiwrite(cast(mask.*cSSE,opts.mapDatatype),'inputReg_SSE',opts.info.map);
                    niftiwrite(cast(mask.*cTstat,opts.mapDatatype),'inputReg_Tstat',opts.info.map);
                else
                    saveImageData(mask.*GLM_Estimate, opts.headers.map, opts.glmlagdir, 'inputReg_beta.nii.gz', datatype);
                    saveImageData(mask.*tmpLag, opts.headers.map, opts.glmlagdir, 'inputReg_lags.nii.gz', datatype);
                    saveImageData(mask.*GLM_lags, opts.headers.map, opts.glmlagdir, 'uncor_inputReg_lags.nii.gz', datatype);
                    saveImageData(mask.*cR2, opts.headers.map, opts.glmlagdir, 'inputReg_R2.nii.gz', datatype);
                    saveImageData(mask.*cSSE, opts.headers.map, opts.glmlagdir, 'inputReg_SSE.nii.gz', datatype);
                    saveImageData(mask.*cTstat, opts.headers.map, opts.glmlagdir, 'inputReg_Tstat.nii.gz', datatype);
                end
                maps.GLM.inputReg_ES = GLM_Estimate;
                maps.GLM.inputReg_lags = tmpLag;
                maps.GLM.uncor_inputReg_lags = GLM_lags;
                maps.GLM.inputReg_R2 = cR2;
                maps.GLM.inputReg_SSE = cSSE;
                maps.GLM.inputReg_Tstat = cTstat;
                
                %if opts.verbose; saveas(gcf,[opts.glmlagdir,'regression_inputReg.fig']); end
        end
        
        %%%%% generate lag-corrected CVR maps %%%%
        if opts.cvr_maps
            disp('Generating lag-adjusted CVR maps based on GLM analysis')
            cCVR = zeros([1 xx*yy*zz]);
            shifted_regr = NaN([length(coordinates) dyn*opts.interp_factor]);
            
            regr = probe;
            for ii = 1:length(coordinates)
                corr_regr = circshift(regr',maxindex(ii));
                if maxindex > 0
                    corr_regr(1,1:maxindex(ii)) = regr(1,1);
                elseif index < 0
                    corr_regr(1,end-maxindex:end) = regr(1,end);
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
                if opts.niiwrite
                    cd(opts.glmCVRdir)
                    niftiwrite(cast(mask.*cCVR,opts.mapDatatype),'optiReg_cCVR',opts.info.map);
                else
                    saveImageData(mask.*cCVR, opts.headers.map, opts.glmCVRdir, 'optiReg_cCVR.nii.gz', datatype);
                end
            else
                if opts.niiwrite
                    cd(opts.glmCVRdir)
                    niftiwrite(cast(mask.*cCVR,opts.mapDatatype),'inputReg_cCVR',opts.info.map);
                else
                    saveImageData(mask.*cCVR, opts.headers.map, opts.glmCVRdir, 'inputReg_cCVR.nii.gz', datatype);
                end
            end
            maps.GLM.CVR.optiReg_cCVR = cCVR;
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
                if opts.niiwrite
                    cd(opts.glmCVRdir)
                    niftiwrite(cast(mask.*cR2,opts.mapDatatype),'optiReg_cR2',opts.info.map);
                    niftiwrite(cast(mask.*cSSE,opts.mapDatatype),'optiReg_cSSE',opts.info.map);
                    niftiwrite(cast(mask.*cTstat,opts.mapDatatype),'optiReg_cTstat',opts.info.map);
                else
                    saveImageData(mask.*cR2, opts.headers.map, opts.glmCVRdir, 'optiReg_cR2.nii.gz', datatype);
                    saveImageData(mask.*cSSE, opts.headers.map, opts.glmCVRdir, 'optiReg_cSSE.nii.gz', datatype);
                    saveImageData(mask.*cTstat, opts.headers.map, opts.glmCVRdir,'optiReg_cTstat.nii.gz', datatype);
                end
                maps.GLM.CVR.optiReg_cR2 = cR2;
                maps.GLM.CVR.optiReg_cSSE = cSSE;
                maps.GLM.CVR.optiReg_cTstat = cTstat;
                
            else
                if opts.niiwrite
                    cd(opts.glmCVRdir)
                    niftiwrite(cast(mask.*cR2,opts.mapDatatype),'inputReg_cR2',opts.info.map);
                    niftiwrite(cast(mask.*cSSE,opts.mapDatatype),'inputReg_cSSE',opts.info.map);
                    niftiwrite(cast(mask.*cTstat,opts.mapDatatype),'inputReg_cTstat',opts.info.map);
                else
                    saveImageData(mask.*cR2, opts.headers.map, opts.glmCVRdir, 'inputReg_cR2.nii.gz', datatype);
                    saveImageData(mask.*cSSE, opts.headers.map, opts.glmCVRdir, 'inputReg_cSSE.nii.gz', datatype);
                    saveImageData(mask.*cTstat, opts.headers.map, opts.glmCVRdir,'inputReg_cTstat.nii.gz', datatype);
                end
            end
            maps.GLM.CVR.inputReg_cR2 = cR2;
            maps.GLM.CVR.inputReg_cSSE = cSSE;
            maps.GLM.CVR.inputReg_cTstat = cTstat;
        end
    end
    disp(['finished running GLM analysis in: ',int2str((cputime-q)/60),' minutes'])
    disp('saving maps in .mat file' )
    maps.newprobe = newprobe;
    save([opts.resultsdir,'lagCVR_maps.mat'], 'maps');
    %   close all
end

