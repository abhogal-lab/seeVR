% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% This implementation was developped by Stefan Rademakers,
% stefan-rademakers@outlook.com
% Sections of this code were contributed by Allen A. Champagne, a.champagne@queensu.ca
% The recursive lag approach is based on the origical code shared by Dr.
% Toshiko Aso:
% https://github.com/aso-toshihiko/BOLDLagMapping_Deperfusioning
% (aso.toshihiko@gmail.com )
%
% <recursiveLag: recursive estimation of perfusion lag structure >
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

function [maps] = recursiveLag(refmask,mask,data,probe,nuisance,opts)
refmask = logical(refmask); mask = logical(mask);
warning('off');
global opts;

% check for statistics toolbox
if license('test','Statistics_toolbox') == 0; opts.pca_analysis = 0; end

maps = struct();
if iscolumn(probe); else; probe = probe'; end
if isempty(nuisance); nuisance = ones(size(probe)); end
test1 = nuisance(1,:); test2 = nuisance(:,1);
if length(test1) > length(test2); nuisance = nuisance'; end; clear test1 test2

% setup default parameters
if isfield(opts,'gpu'); else; opts.gpu = 0; end                           %use gpu
if isfield(opts,'niiwrite'); else; opts.niiwrite = 0; end                 % depending on how data is loaded this can be set to 1 to use native load/save functions
if isfield(opts,'plot'); else; opts.plot = 0; end                         % show or hide selected plots
if isfield(opts,'prewhite'); else; opts.prewhite = 0; end                 % zero mean and unit variance of data. Normally leads to bad results
if isfield(opts,'interp_factor'); else; opts.interp_factor = 4; end       % factor by which to temporally interpolate data. Better for picking up lags between TR
if isfield(opts,'load_probe'); else; opts.load_probe = 0; end             % saves time for creating regressor if one with the correct length already exists
if isfield(opts,'glm_model'); else; opts.glm_model = 0; end               % perform GLM based lag regression. This works well for high temporal resolution data
if isfield(opts,'uni'); else; opts.uni = 0; end                           % if this is set to 1 negative correlations will be ignored
if isfield(opts,'useGPU'); else; opts.useGPU = 0; end                     % alot of cards may not have enough memory for this so default is no


% setup main default parameters
if isfield(opts,'rescale_probe'); else; opts.rescale_probe = 1; end       % rescaling may be helful for initial refinement run
if isfield(opts,'trace_corr'); else; opts.trace_corr = 1; end             % perform additional correlation with gas trace (on top of optimized regressor)
if isfield(opts,'refine_regressor'); else; opts.refine_regressor = 1; end % refine BOLD regressor. If '0' a straight correlation will be done with probe
if isfield(opts,'pca_analysis'); else; opts.pca_analysis = 1; end         % PCA analysis to optimize BOLD regressor
if isfield(opts,'corr_model'); else; opts.corr_model = 1; end             % perform correlation based analysis
if isfield(opts,'norm_regr'); else; opts.norm_regr = 0; end               % normalize the regressor for correlation analysis between 0-1

% hyperparameters for recursive lag mapping
if isfield(opts,'lim'); else; opts.lim = 3; end                           % range of locally calculated cross-correlations
if isfield(opts,'lim_s'); else; opts.lim_s = 1; end                       % range of indexes whose corresponding local cross-correlations are saved
if isfield(opts,'THR'); else; opts.THR = 0.2; end                         % cross-correlation acceptance threshold
if isfield(opts,'comp'); else; opts.comp = 1; end                         % enabling a comparison between recursive and 'classic' lag mapping if set to 1

if opts.niiwrite
    if isfield(opts.info,'rts'); else; opts.info.rts = opts.info.ts; end
end

% important parameters for results
if isfield(opts,'corr_thresh'); else; opts.corr_thresh = 0.7; end         % threshold by which to accept correlated voxels during the refinement of the BOLD regressor

% refinement
if isfield(opts,'lowerlagthresh'); else; opts.lowerlagthresh = -2; end    % lower threshold for refinement (generaly -2 to 0)
if isfield(opts,'upperlagthresh'); else; opts.upperlagthresh = 2; end     % upper threshold for refinement (generaly 0 to +2)

% account for interpolation factor
opts.adjlowerthresh = opts.lowerlagthresh*opts.interp_factor;
opts.adjupperthresh = opts.upperlagthresh*opts.interp_factor;

% correlation
if isfield(opts,'lowlag'); else; opts.lowlag = -3; end                    % lower threshold for correlation (generaly -3 to 0)
if isfield(opts,'highlag'); else; opts.highlag = 40; end                  % upper threshold for correlation (in healthy up to 20-40, in disease 60-90)

% account for interpolation factor
opts.adjlowlag = opts.lowlag*opts.interp_factor;                          % setups lower lag limit; negative for misalignment and noisy correlation
opts.adjhighlag = opts.highlag*opts.interp_factor;                        % setups upper lag limit; allow for long lags associated with pathology

%setup save directories
    if isfield(opts,'resultsdir'); else; opts.resultsdir = pwd; end
    opts.corrlagdir = fullfile(opts.resultsdir,'corrLAG'); mkdir(opts.corrlagdir);
    
%optional image outputs
cd(opts.resultsdir);
datatype = 16;
[xx, yy, zz, dyn] = size(data);
orig_regr = probe;

t = cputime;

if opts.useGPU
if ~(parallel.gpu.GPUDevice.isAvailable)
    fprintf(['\n\t**GPU not available. Stopping.**\n']);
    opts.GPUcorr = 0;
    return;
else
    dev = gpuDevice;
    fprintf(...
    'GPU detected (%s, %d multiprocessors, Compute Capability %s)',...
    dev.Name, dev.MultiprocessorCount, dev.ComputeCapability);
opts.GPUcorr = 1;
end
end


%% grab coordinates
[orig_voxel_ts, coordinates] = grabTimeseries(data, mask);    % WB coordinates
[gm_voxel_ts, gmcoordinates] = grabTimeseries(data, refmask); % GM coordinates 

% perform pre-whitening
if opts.prewhite
    gm_voxel_ts = gm_voxel_ts'; parfor ii=1:size(gmcoordinates,1); [gm_voxel_ts(:,ii), ~, ~] = prewhiten(gm_voxel_ts(:,ii)); end; gm_voxel_ts = gm_voxel_ts';
    pw_voxel_ts = orig_voxel_ts'; parfor ii=1:size(coordinates,1); [pw_voxel_ts(:,ii), ~, ~] = prewhiten(pw_voxel_ts(:,ii)); end; pw_voxel_ts = pw_voxel_ts';
    [probe, ~, ~] = prewhiten(probe);
end

% interpolate variables
probe = interp(probe,opts.interp_factor); % interpolate data by a factor of the number of slices for more accurate timing
orig_regr = interp(orig_regr,opts.interp_factor); % save input regressor to generate corrected CVR maps
for ii=1:size(nuisance,2); np(:,ii) = demeanData(rescale(interp(nuisance(:,ii),opts.interp_factor))); end
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

% grabbing WB timeseries
tic
wb_voxel_ts=zeros([length(coordinates) opts.interp_factor*dyn]);
clean_voxel_ts=zeros([length(coordinates) opts.interp_factor*dyn]);

% For lags you want to use wb_voxel_ts (accounting for prewhitening). 
% For CVR maps always regress original regressor with clean timeseries
if opts.prewhite
    parfor ii = 1:length(coordinates)
        clean_voxel_ts(ii,:) = interp(orig_voxel_ts(ii,:),opts.interp_factor);
        wb_voxel_ts(ii,:) = interp(pw_voxel_ts(ii,:),opts.interp_factor);
    end
else
    parfor ii = 1:length(coordinates)
        clean_voxel_ts(ii,:) = interp(orig_voxel_ts(ii,:),opts.interp_factor);
    end
    wb_voxel_ts = clean_voxel_ts;
end

%% Recursive method
disp('performing recursive lag mapping')
% define the parameters for the whole set of possible lags
lags = opts.adjlowlag:opts.adjhighlag;            % 'global' range of possible lags
gmid = find(lags==0);                             % the index of lag 0

% define the parameters for the local set of lags calculated per iteration
lmid = opts.lim+1;                        % middle of the limited range of calculated lags
lrange = -opts.lim_s:opts.lim_s;          % limited range of indexes whose lag is saved
lrange_s = lmid+lrange;                   % limited saved range shifted to middle of limited calculated range

% initialization
r_map = zeros([xx*yy*zz length(lags)]);
lag_map = NaN([xx*yy*zz length(lags)]);

% WB cross-correlation to the probe across the limited range of lags
a2=mat2cell(wb_voxel_ts,ones(1,size(wb_voxel_ts,1)),size(wb_voxel_ts,2));
b2=cellfun(@(a2) xcorr(a2,newprobe,opts.lim,'coeff'),a2,'uniformoutput',false);
corr_probe = cell2mat(b2); corr_probe = corr_probe';
if opts.uni
    [R, index] = max(corr_probe); % negative correlations will not weigh more than positive ones
else
    [R, index] = max(abs(corr_probe)); % absolute max to take care of negative correlation
end
R(R<opts.THR) = 0; % remove low r-values

% comparison between the adj Aso method and a simple cross correlation to the initial regressor
if opts.comp
    check_b2=cellfun(@(a2) xcorr(a2,newprobe,'coeff'),a2,'uniformoutput',false);
    check_corr_probe = cell2mat(check_b2); check_corr_probe = check_corr_probe';
    idx = find(lags<opts.adjlowerthresh | lags>opts.adjupperthresh);
    check_corr_probe(idx,:)=[]; clear idx
    if opts.uni
        [check_r, ~] = max(check_corr_probe); % negative correlations will not weigh more than positive ones
    else
        [check_r, ~] = max(abs(check_corr_probe)); % absolute max to take care of negative correlation
    end
end

% perform first iteration to allocate lag 0 (and other lags within lrange_s) and allocate its respective r-values
r_map(coordinates(ismember(index,lrange_s)),gmid) = R(ismember(index,lrange_s));
for ind = lrange_s
    add = lrange(lrange_s==ind);
    lag_map(coordinates(ismember(index,ind)),gmid) = lags(gmid+add);
end

% perform consecutive iterations, both in positive and negative lags
Uregr = zeros(size(newprobe)); Lregr = zeros(size(newprobe)); % initialization

%% positive p
it = 1; 
for p = 1:opts.adjhighlag-1
    disp('it: '+string(it)+'/'+string(length(opts.adjlowlag:opts.adjhighlag)-1))
    it = it + 1;
    
    % use original index if iteration is too low to update regressor
    if p<=opts.lim
        Uindex = index;
        Lindex = index;
    end
    
    % save temporal data of WB voxels
    data = zeros(size(wb_voxel_ts));
    data(:,1+opts.lim+p:end-opts.lim-p) = wb_voxel_ts(:,1+opts.lim+p:end-opts.lim-p);
    a2=mat2cell(data,ones(1,size(data,1)),size(data,2));
    
    % form a new shifted regressor of the next positive lag step
    if opts.rescale_probe;
    regr = rescale(mean(wb_voxel_ts(ismember(Uindex,lmid+1),:),1));
    else
    regr = mean(wb_voxel_ts(ismember(Uindex,lmid+1),:),1);
    end
    figure(20); hold on; plot(regr);
    % shift the regressor such that its temporal position corresponds to the BOLD data
    Uregr(1+p:end) = regr(1:end-p); plot(rescale(Uregr));
    
    % find maximum cross-correlation across the limited range of calculated lags
    b2=cellfun(@(a2) xcorr(a2,Uregr,opts.lim,'coeff'),a2,'uniformoutput',false);
    corr_probe = cell2mat(b2); corr_probe = corr_probe';
    if opts.uni
        [R, Uindex] = max(corr_probe); % negative correlations will not weigh more than positive ones
    else
        [R, Uindex] = max(abs(corr_probe)); % absolute max to take care of negative correlation
    end
    %R(R<opts.THR) = 0; % remove low r-values
    
    % save lags and r that are within lrange_s
    r_map(coordinates(ismember(Uindex,lrange_s)),gmid+p) = R(ismember(Uindex,lrange_s));
    for ind = lrange_s
        add = lrange(lrange_s==ind);
        lag_map(coordinates(ismember(Uindex,ind)),gmid+p) = lags(gmid+add+p);
    end

    % negative p
    if p<=abs(opts.adjlowlag)-1
        disp('it: '+string(it)+'/'+string(length(opts.adjlowlag:opts.adjhighlag)-1))
        it = it + 1;
    
        % form a new shifted regressor of the next negative lag step
        regr = mean(wb_voxel_ts(ismember(Lindex,lmid-1),:),1);
        
        % shift the regressor such that its temporal position corresponds to the BOLD data
        Lregr(1:end-p) = regr(1+p:end); plot(rescale(Lregr));
        
        % find maximum cross-correlation across the limited range
        b2=cellfun(@(a2) xcorr(a2,Lregr,opts.lim,'coeff'),a2,'uniformoutput',false);
        corr_probe = cell2mat(b2); corr_probe = corr_probe';
        if opts.uni
            [R, Lindex] = max(corr_probe); % negative correlations will not weigh more than positive ones
        else
            [R, Lindex] = max(abs(corr_probe)); % absolute max to take care of negative correlation
        end
        R(R<opts.THR) = 0; % remove low r-values
        
        % save lags and r that are within lrange_s
        r_map(coordinates(ismember(Lindex,lrange_s)),gmid-p) = R(ismember(Lindex,lrange_s));
        for ind = lrange_s
            add = lrange(lrange_s==ind);
            lag_map(coordinates(ismember(Lindex,ind)),gmid-p) = lags(gmid+add-p);
        end
    end
end

% assign only the lag and r values of the highest r
[max_r_map, max_idx] = max(r_map,[],2);
lag_map_max = NaN([xx*yy*zz 1]);
for ii = 1:length(coordinates)
    lag_map_max(coordinates(ii)) = lag_map(coordinates(ii),max_idx(coordinates(ii)));
end
lag_map = lag_map_max; r_map = max_r_map;
tmpLag = opts.TR*(lag_map/opts.interp_factor); % set Lags back to seconds

% calculate and display success (and improvement if opts.comp=1)
succesfull = sum(~isnan(lag_map(coordinates)));
perc = 100*succesfull/length(coordinates);
disp(string(perc)+'% of voxels are allocated to lag_map')
if opts.comp
    perc = 100*sum(max_r_map(coordinates)>=check_r')/sum(check_r>0);
    disp(string(perc)+'% of voxels have an improved r-value w.r.t. the original cross-correlation')
end

% reshape and save lag and r maps
r_map = reshape(r_map,[xx yy zz]);
tmpLag = reshape(tmpLag,[xx yy zz]);
maps.XCORR.recursiveLag_map = mask.*tmpLag;
maps.XCORR.recursiveR_map = mask.*r_map;
if opts.niiwrite
    niftiwrite(cast(mask.*tmpLag,opts.mapDatatype),'recursiveLag',opts.info.map);
    niftiwrite(cast(mask.*r_map,opts.mapDatatype),'recursiveR_map',opts.info.map);
else
    saveImageData(mask.*tmpLag, opts.headers.map, opts.corrlagdir, 'recursiveLag.nii.gz', datatype);
    saveImageData(mask.*r_map, opts.headers.map, opts.corrlagdir, 'recursiveR_map.nii.gz', datatype);
end
end

