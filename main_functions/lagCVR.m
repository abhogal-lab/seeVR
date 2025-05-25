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
try
    warning('off');
    warning('off', 'MATLAB:rankDeficientMatrix');

    spmd
        warning('off', 'MATLAB:rankDeficientMatrix');
    end
catch
end
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
    opts.motioncorr = 0.3;
    [np, ~] = prepNuisance(nuisance,probe, opts);
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
if isfield(opts,'trace_corr'); else; opts.trace_corr = 1; end             %perform additional correlation with gas trace (on top of optimized regressor)
if isfield(opts,'refine_regressor'); else; opts.refine_regressor = 1; end %refine BOLD regressor. If '0' a straight correlation will be done with probe
if isfield(opts,'corr_model'); else; opts.corr_model = 1; end             %perform correlation based analysis
if isfield(opts,'cvr_maps'); else; opts.cvr_maps = 1; end                 %generate CVR maps based on regression of BOLD signal with gas probe
if isfield(opts,'eff_probe'); else; opts.eff_probe = 0; end               %uses effective probe to calculate probe response
if isfield(opts,'glm_model'); else; opts.glm_model = 0; end               %perform GLM based lag regression. This works well for high temporal resolution data
if isfield(opts,'uni'); else; opts.uni = 0; end                           %if this is set to 1 negative correlations will be ignored
if isfield(opts,'norm_regr'); else; opts.norm_regr = 0; end               %normalize the regressor for correlation analysis between 0-1
if isfield(opts,'robust'); else; opts.robust = 0; end                     %calculate robust lag and CVR
if isfield(opts,'refine_lag'); else; opts.refine_lag = 1; end             %When set to 1, lag calculation will reprocess voxels with clipped values by creating a new mean time series that averages data from neighboring voxels
if isfield(opts,'win_size'); else; opts.win_size = 1; end                 %the number of voxels to consider around the voxel of interest when opts.refine_lag = 1;
if isfield(opts,'passes'); else; opts.passes = 10; end                    %how many iterations for lag_map refinement
if isfield(opts,'medfilt_maps'); else; opts.medfilt_maps = 1; end         %cleans up lag map with median filter but does not change statistical maps. May improve lag-corrected CVR
if isfield(opts,'filloutliers'); else; opts.filloutliers = 1; end         %spike removal
if isfield(opts,'figdir'); else; opts.figdir = opts.resultsdir; end        
if opts.niiwrite
    if isfield(opts.info,'rts'); else; opts.info.rts = opts.info.ts; end
end

%important parameters for results
if isfield(opts,'corr_thresh'); else; opts.corr_thresh = 0.7; end         %threshold by which to accept correlated voxels during the refinement of the BOLD regressor
%refinement
if isfield(opts,'lowerlagthresh'); else; opts.lowerlagthresh = -2; end    %lower threshold for refinement (generaly -2 to 0)
if isfield(opts,'upperlagthresh'); else; opts.upperlagthresh = 2; end     %upper threshold for refinement (generaly 0 to +2)

%account for interpolation factor
opts.lowerthresh = opts.lowerlagthresh*opts.interp_factor;
opts.upperthresh = opts.upperlagthresh*opts.interp_factor;

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
if opts.filloutliers
    disp('removing outliers...')
    [orig_voxel_ts, coordinates] = grabTimeseries(data, mask);
    B = filloutliers(orig_voxel_ts', 'spline', 'mean');
    B = B';
    temp = zeros([numel(mask) size(data,4)]);
    temp(coordinates,:) = B; 
    data = reshape(temp, size(data)); clear temp; clear B;
else
    % WB coordinates
    [orig_voxel_ts, coordinates] = grabTimeseries(data, mask);
end

% GM coordinates
[gm_voxel_ts, gmcoordinates] = grabTimeseries(data, refmask);

if opts.prewhite
    gm_voxel_ts = gm_voxel_ts'; parfor ii=1:size(gmcoordinates,1); [gm_voxel_ts(:,ii), ~, ~] = prewhiten(gm_voxel_ts(:,ii)); end; gm_voxel_ts = gm_voxel_ts';
    pw_voxel_ts = orig_voxel_ts'; parfor ii=1:size(coordinates,1); [pw_voxel_ts(:,ii), ~, ~] = prewhiten(pw_voxel_ts(:,ii)); end; pw_voxel_ts = pw_voxel_ts';
    [probe, ~, ~] = prewhiten(probe);
end


%% to save time, check if there is an existing probe
if opts.refine_regressor
    probename = 'final_probe.mat'
    if exist(probename) && opts.load_probe
        newprobe = load(probename); newprobe = newprobe.newprobe;
        disp('found existing probe; skipping generation of BOLD regressor')
    else
        [newprobe] = optimizeRegressor(probe, data, refmask, opts)
        save(fullfile(opts.resultsdir,'final_probe.mat'), 'newprobe');
        disp('Finished creating optimized regressor')
        clear a GMmask yes score coeff gm_voxel_ts 
    end
else
    newprobe = probe;
    opts.trace_corr = 0;
    [~,lags] = xcorr(newprobe,wb_voxel_ts(1,:),'coeff');
    idx = find(lags<=opts.adjlowlag | lags>=opts.adjlowlag);
    lags(idx)=[];
end
%% for final round, take the refined regressor and do final
% correlations between it and the original signal to estimate the peak time delay,
% correlation for the whole brain

%interpolate variables
newprobe = interp(newprobe,opts.interp_factor); %interpolate data by a factor of the number of slices for more accurate timing
orig_regr = interp(orig_regr,opts.interp_factor); %save input regressor to generate corrected CVR maps
probe = orig_regr;

%%% grabbing WB timeseries
tic
wb_voxel_ts=zeros([length(coordinates) opts.interp_factor*dyn]);

%For lags you want to use wb_voxel_ts (accounting for prewhitening. For CVR
%maps always regress original regressor with clean timeseries

if opts.cvr_maps && opts.prewhite
    opts.prewhite = 0;
    disp('prewhitening is not compatible with CVR mapping')
    disp('CVR option takes priotity... for lag mapping using prewhitening, set opts.cvr_maps = 1 + opts.prewhite = 1')
end

if opts.prewhite
    parfor ii = 1:length(coordinates)
        wb_voxel_ts(ii,:) = interp(pw_voxel_ts(ii,:),opts.interp_factor);
    end
else
    parfor ii = 1:length(coordinates)
        wb_voxel_ts(ii,:) = interp(orig_voxel_ts(ii,:),opts.interp_factor);
    end
end; clear ip


if opts.corr_model

    if opts.trace_corr; qq = [1 2]; else; qq = 1; end %correlate probe y/n

    switch numel(qq)
        case 1
            index_map = zeros([1 numel(mask)]);  %setup the index mask to properly keep track of replaced lag voxels
            rvec = zeros([1, length(coordinates)]);
            index = zeros([2 size(wb_voxel_ts,1)]);
        case 2
            index_map = zeros([2 numel(mask)]);  %setup the index mask to properly keep track of replaced lag voxels
            rvec = zeros([2, length(coordinates)]);
            index = zeros([2 size(wb_voxel_ts,1)]);
    end

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
        corr_probe=[];a2=[];
        %matrix cross correlation with probe
        a2=mat2cell(wb_voxel_ts,ones(1,size(wb_voxel_ts,1)),size(wb_voxel_ts,2)); %turn into cells so treat all rows independently

        b2=cellfun(@(a2) xcorr(a2,regr,opts.adjhighlag,'none'),a2,'uniformoutput',false);
        corr_probe=cell2mat(b2); %reformat into correlation matrix
        corr_probe = corr_probe';
        %use highlag to restrict the number of correlations that need to be done
        [~,lags] = xcorr(wb_voxel_ts(1,:),regr,opts.adjhighlag,'none');  %%%%%%

        %save correlations over time only if this map is to be saved (saves

        if opts.save_rts
            corr_ts = zeros([xx*yy*zz size(corr_probe,1)]);
            corr_ts(coordinates,:) = corr_probe';
            corr_ts(isnan(corr_ts)) = 0;
            corr_ts = reshape(corr_ts, [xx yy zz size(corr_probe,1)]);
        end

        % trim correlation probe
        idx = find(lags<=opts.adjlowlag | lags>=opts.adjhighlag);
        %remove lag values lower that what is possible
        lags(idx)=[];
        corr_probe(idx,:)=[];

        if opts.uni %add option to pick closest probe to zero-lag
            [~, index_map(pp,coordinates)] = max(corr_probe); %negative correlations will not weigh more than positive ones
            [~, index(pp,:)] = max(corr_probe); %negative correlations will not weigh more than positive ones
        else
            [~, index_map(pp,coordinates)] = max(abs(corr_probe)); %Absolute max to take care of negative correlation
            [~, index(pp,:)] = max(abs(corr_probe)); %Absolute max to take care of negative correlation
        end
        %fill correlation map
        parfor ii = 1:length(coordinates)
            rvec(pp,ii) = corr_probe(index(pp,ii),ii);
        end

        r_map = zeros([xx*yy*zz 1]); lag_map = zeros([1 xx*yy*zz]);
        r_map(coordinates) = rvec(pp,:);
        r_map = reshape(r_map,[xx yy zz]);

        %fill lag map
        lag_map(1,coordinates) = lags(1,index_map(pp,coordinates));
        %        tmpLag = opts.TR*(lag_map/opts.interp_factor); tmpLag(1,coordinates) = tmpLag(1,coordinates) + abs(min(tmpLag(:))); %puts lag maps back into seconds
        lag_map = reshape(lag_map,[xx yy zz]);
        %       tmpLag = reshape(tmpLag,[xx yy zz]);

        %identify clipped lag values and use neighboring information to
        %improve correlation
        if opts.refine_lag
            iter = 1;
            passes = 0;
            maxlag = max(lag_map(:));

            while iter
                passes = passes + 1
                max_map = zeros(size(lag_map));
                max_map(lag_map == maxlag) = 1;

                %get new coordinates
                [~, newcoordinates] = grabTimeseries(data, max_map);

                i = find(max_map); % find nonzero values in M
                [X,Y,Z] = ind2sub(size(max_map), i); clear i
                newTS = zeros([length(X) length(regr)]);

                for kk=1:length(X)
                    try
                        %extract timeseries and neighboring timeseries
                        X_rng = X(kk)-opts.win_size:X(kk)+opts.win_size;
                        Y_rng = Y(kk)-opts.win_size:Y(kk)+opts.win_size;
                        Z_rng = Z(kk)-opts.win_size:Z(kk)+opts.win_size;
                        %remove possiblity to have slices outside FOV
                        X_rng(X_rng > size(data,1)) = []; X_rng(X_rng < 1) = [];
                        Y_rng(Y_rng > size(data,2)) = []; Y_rng(Y_rng < 1) = [];
                        Z_rng(Z_rng > size(data,3)) = []; Z_rng(Z_rng < 1) = [];

                        %extract timeseries
                        tmp =  data(X_rng,Y_rng,Z_rng,:);
                        tmp = reshape(tmp, [size(tmp,1)*size(tmp,2)*size(tmp,3), size(tmp,4)]);
                        tmp(find(tmp(:,1) == 0),:) = [];
                        data(X(kk),Y(kk),Z(kk),:) = mean(tmp); %update data matrix for next iteration
                        newTS(kk,:) = interp(mean(tmp),opts.interp_factor);
                    catch
                        disp('error; skipping voxel')
                    end
                end

                clear tmp
                corr_probe2=[];a2=[];b2=[];

                %matrix cross correlation with probe
                a2=mat2cell(newTS,ones(1,size(newTS,1)),size(newTS,2)); %turn into cells so treat all rows independently
                b2=cellfun(@(a2) xcorr(a2,regr,opts.adjhighlag,'none'),a2,'uniformoutput',false);
                corr_probe=cell2mat(b2); %reformat into correlation matrix
                corr_probe = corr_probe';
                corr_probe(idx,:)=[];

                %create lag maps
                rvec2 = zeros([1, length(X)]);
                index2 = zeros([1 size(newTS,1)]);

                if opts.uni
                    [~, index2] = max(corr_probe); %negative correlations will not weigh more than positive ones
                else
                    [~, index2] = max(abs(corr_probe)); %Absolute max to take care of negative correlation
                end

                parfor ii = 1:length(X)
                    rvec2(1,ii) = corr_probe(index2(1,ii),ii);
                end

                lag_map = lag_map(:);
                r_map = r_map(:);
                r_map(newcoordinates) = rvec2;
                r_map = reshape(r_map,[xx yy zz]);
                lag_map(newcoordinates) = lags(1,index2);
                lag_map = reshape(lag_map,[xx yy zz]);

                %update the index for the corrected map
                index_map(pp,newcoordinates) = index2;
                perc = 100*(length(index2)/length(coordinates))
                disp([int2str(perc), ' percent of voxels have clipped lag values'])
                if perc > 2
                    if passes < opts.passes
                        continue;
                    else; iter = 0;
                        disp('exceeded the maximum allowed passes')
                        disp('... to increase passes set opts.passes to a higher value')
                    end
                else
                    iter = 0;
                end

                clear rvec2 index2;
            end
        end

        %save lag and r maps
        tmpLag = opts.TR*(lag_map/opts.interp_factor).*mask;
        if min(tmpLag(:)) < 0; tmpLag = tmpLag + abs(min(tmpLag(:))); end

        if opts.medfilt_maps
            tmpLag = medfilt3(tmpLag);
            lag_map = medfilt3(lag_map);
        end
        
        savedir = opts.corrlagdir;
        switch pp
            case 1

                saveMap(cast(mask.*tmpLag,opts.mapDatatype), savedir, 'hemodynamic_lag_map_refined_probe', opts.info.map, opts);
                saveMap(cast(mask.*lag_map,opts.mapDatatype), savedir, 'raw_hemodynamic_lag_map_refined_probe', opts.info.map, opts);
                saveMap(cast(mask.*r_map,opts.mapDatatype), savedir, 'r_map_refined_probe', opts.info.map, opts);
                
                maps.XCORR.lag_opti = tmpLag;
                maps.XCORR.uncorrlag_opti = lag_map;
                maps.XCORR.r_opti = r_map;

                if opts.save_rts
                    if opts.niiwrite
                        cd(opts.corrlagdir);
                        disp('saving correlation timeseries based on optimized regressor')
                        niftiwrite(cast(corr_ts,opts.mapDatatype),'correlation_timeseries_refined_probe',opts.info.rts);
                        clear Trcorr_ts;
                    else
                        disp('saving correlation timeseries based on optimized regressor')
                        saveImageData(corr_ts, headers.ts, opts.corrlagdir,  'correlation_timeseries_refined_probe.nii.gz', datatype)
                        clear Trcorr_ts;
                    end
                end
                if opts.robust
                    LAG(1,:,:,:) = mask.*tmpLag;
                    RL(1,:,:,:) = mask.*r_map;
                end
                clear tmpLag lag_map r_map;
            case 2

                saveMap(cast(mask.*tmpLag,opts.mapDatatype), savedir, 'hemodynamic_lag_map_input_probe', opts.info.map, opts);
                saveMap(cast(mask.*lag_map,opts.mapDatatype), savedir, 'raw_hemodynamic_lag_map_input_probe', opts.info.map, opts);
                saveMap(cast(mask.*r_map,opts.mapDatatype), savedir, 'r_map_input_probe', opts.info.map, opts);
               
                maps.XCORR.lag_input = tmpLag;
                maps.XCORR.uncorrlag_input = lag_map;
                maps.XCORR.r_input = r_map;

                if opts.save_rts
                        saveMap(cast(corr_ts,opts.mapDatatype), savedir, 'correlation_timeseries', opts.info.map, opts);
                        disp('saving correlation timeseries')
                        clear Trcorr_ts;
                end
                if opts.robust
                    LAG(2,:,:,:) = mask.*mask.*tmpLag;
                    RL(2,:,:,:) = mask.*r_map;
                end
                clear tmpLag lag_map r_map;
        end
        clear Trcorr_ts orig_voxel_ts newTS a2 b2 corr_probe;
    end

    index = index_map(:, coordinates);
    disp(['Lag, regressor and r_maps were created in: ',int2str((cputime-t)/60),' minutes'])

    clear m corr_ts regr index_map;
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
        niftiwrite(cast(logical(mask).*robustIR,opts.mapDatatype),'robust_hemodymic_lag_map_r',opts.info.map);
    else
        saveImageData(mask.*robustIR, opts.headers.map, opts.corrlagdir,'robust_hemodymic_lag_map_r.nii.gz', datatype);
    end
    maps.XCORR.robust_lag_r = robustIR;
    clear LAG RL
else
    clear LAG RL
end


%% Generate lag adjusted CVR maps

if opts.cvr_maps
    %% normalize probe and refined regressor, determine effective CO2 trace
    if opts.trace_corr; index(2,:) = []; end % if 2 runs are done (input and optimized probe) ignore the input indices since optimization is generally better
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
                saveas(gcf,fullfile(opts.figdir,'effective_probes.fig'));
                disp('Generating base maps using effective probe trace')
                regr = eff_probe;
        end

        bCVR = zeros([1 xx*yy*zz]); bR2 = zeros([1 xx*yy*zz]); bSSE = zeros([1 xx*yy*zz]); bTstat = zeros([1 xx*yy*zz]);

        %create base CVR map
        A = regr;
        C = [ones([length(A) 1]) A]; clear A
        regr_coef = C\wb_voxel_ts';

        bCVR(1,coordinates) = regr_coef(2,:); %extract slope
        bCVR(bCVR > 30) = 0; bCVR(bCVR < -30) = 0; %cleanup base CVR map
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
                    niftiwrite(cast(mask.*bCVR,opts.mapDatatype),'lin_regr_CVR_map',opts.info.map);
                    niftiwrite(cast(mask.*bR2,opts.mapDatatype),'lin_regr_CVR_R2_map',opts.info.map);
                    niftiwrite(cast(mask.*bTstat,opts.mapDatatype),'lin_regr_CVR_Tstat_map',opts.info.map);
                else
                    saveImageData(mask.*bCVR, opts.headers.map, opts.corrCVRdir, 'lin_regr_CVR_map.nii.gz', datatype);
                    saveImageData(mask.*bR2, opts.headers.map, opts.corrCVRdir,'lin_regr_CVR_R2_map.nii.gz', datatype);
                    saveImageData(mask.*bTstat, opts.headers.map, opts.corrCVRdir,'lin_regr_CVR_Tstat_map.nii.gz', datatype);
                end
                maps.XCORR.CVR.bCVR = bCVR;
                maps.XCORR.CVR.bR2 = bR2;
                maps.XCORR.CVR.bSSE = bSSE;
                maps.XCORR.CVR.bTstat = bSSE;

                disp('Base CVR, r2, Tstat and SSE maps were created using endtidal regressor')
                CVR(1,:,:,:) = mask.*bCVR;
                TSTAT(1,:,:,:) = mask.*bTstat;
                RC(1,:,:,:) = mask.*bR2;
                clear bCVR bR2 bT bTstat bSSE SSE SST R2 X Y r regr_coef A C s

            case 2 %effective probe
                if opts.niiwrite
                    opts.corrCVRdir
                    niftiwrite(cast(mask.*bCVR,opts.mapDatatype),'lin_regr_CVR_effective_probe_map',opts.info.map);
                    niftiwrite(cast(mask.*bR2,opts.mapDatatype),'lin_regr_CVR_effective_probe_R2_map',opts.info.map);
                    niftiwrite(cast(mask.*bTstat,opts.mapDatatype),'lin_regr_CVR_effective_probe_Tstat_map',opts.info.map);
                else
                    saveImageData(mask.*bCVR, opts.headers.map, opts.corrCVRdir, 'lin_regr_CVR_effective_probe_map.nii.gz', datatype);
                    saveImageData(mask.*bR2, opts.headers.map, opts.corrCVRdir,'lin_regr_CVR_effective_probe_R2_map.nii.gz', datatype);
                    saveImageData(mask.*bTstat, opts.headers.map, opts.corrCVRdir,'lin_regr_CVR_effective_probe_Tstat_map.nii.gz', datatype);
                end
                maps.XCORR.CVR.bCVR_eff = bCVR;
                maps.XCORR.CVR.bR2_eff = bR2;
                maps.XCORR.CVR.bSSE_eff = bSSE;
                maps.XCORR.CVR.bTstat_eff = bSSE;

                disp('Base CVR, r2, Tstat and SSE maps were created using effective endtidal regressor')
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
            if index(ii) > 0
                corr_regr(1,1:index(ii)) = regr(1,1); %formerly NaN
            elseif index(ii) < 0
                corr_regr(1,end-index(ii):end) = regr(1,end); %formerly NaN
            else
                %do nothing because index is zero = zero shift
            end
            shifted_regr(ii,:) = corr_regr;
        end

        regr_coef = zeros([2 length(coordinates)]);
        parfor ii=1:length(coordinates)
            A = shifted_regr(ii,:);
            B = wb_voxel_ts(ii,:); B(isnan(A)) = []; A(isnan(A)) = [];
            C = [ones([length(A) 1]) A'];
            regr_coef(:,ii)= C\B';
        end

        cCVR(1,coordinates) = regr_coef(2,:); %extract slope
        cCVR(cCVR > 30) = 0; cCVR(cCVR < -30) = 0; %cleanup base CVR map
        cCVR = reshape(cCVR, [xx yy zz]);

        %calculate statistics

        SSE = zeros([1 size(wb_voxel_ts,1)]); SST = zeros([1 size(wb_voxel_ts,1)]); cT = zeros([1 size(wb_voxel_ts,1)]);
        parfor ii=1:size(wb_voxel_ts,1)
            A = shifted_regr(ii,:);
            X = wb_voxel_ts(ii,:); X(isnan(A)) = []; A(isnan(A)) = [];
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
            case 1 %endtidal
                if opts.niiwrite
                    cd(opts.corrCVRdir);
                    niftiwrite(cast(mask.*cCVR,opts.mapDatatype),'lag_corrected_CVR_map',opts.info.map);
                    niftiwrite(cast(mask.*cR2,opts.mapDatatype),'lag_corrected_CVR_R2_map',opts.info.map);
                    niftiwrite(cast(mask.*cTstat,opts.mapDatatype),'lag_corrected_CVR_Tstat_map',opts.info.map);
                else
                    saveImageData(mask.*cCVR, opts.headers.map, opts.corrCVRdir,'lag_corrected_CVR_map.nii.gz', datatype);
                    saveImageData(mask.*cR2, opts.headers.map, opts.corrCVRdir,'lag_corrected_CVR_R2_map.nii.gz', datatype);
                    saveImageData(mask.*cTstat, opts.headers.map, opts.corrCVRdir,'lag_corrected_CVR_Tstat_map.nii.gz', datatype);
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
                        niftiwrite(cast(lagregressor,opts.mapDatatype),'correlation_timeseries_lagregressor_map',opts.info.rts);
                    else
                        saveImageData(lagregressor, opts.headers.ts, opts.corrCVRdir, 'correlation_timeseries_lagregressor_map.nii.gz', datatype);
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
                    niftiwrite(cast(mask.*cCVR,opts.mapDatatype),'lag_corrected_effective_CVR_map',opts.info.map);
                    niftiwrite(cast(mask.*cR2,opts.mapDatatype),'lag_corrected_effective_CVR_R2_map',opts.info.map);
                    niftiwrite(cast(mask.*cTstat,opts.mapDatatype),'lag_corrected_effective_CVR_Tstat_map',opts.info.map);
                else
                    saveImageData(mask.*cCVR, opts.headers.map, opts.corrCVRdir,'lag_corrected_effective_CVR_map.nii.gz', datatype);
                    saveImageData(mask.*cR2, opts.headers.map, opts.corrCVRdir,'lag_corrected_effective_CVR_R2_map.nii.gz', datatype);
                    saveImageData(mask.*cTstat, opts.headers.map, opts.corrCVRdir,'lag_corrected_effective_CVR_Tstat_map.nii.gz', datatype);
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
                        niftiwrite(cast(lagregressor,opts.mapDatatype),'correlation_timeseries_effective_lagregressor_map',opts.info.rts);
                    else
                        saveImageData(lagregressor, opts.headers.ts, opts.corrCVRdir, 'correlation_timeseries_effective_lagregressor_map.nii.gz', datatype);
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
            niftiwrite(cast(mask.*robustIT,opts.mapDatatype),'CVR_based_on_highest_Tstat_map',opts.info.map);
        else
            saveImageData(mask.*robustIT, opts.headers.map, opts.corrCVRdir,'CVR_based_on_highest_Tstat_map.nii.gz', datatype);
        end
        maps.XCORR.CVR.robustCVR_TSTAT = robustIT;
    end
    if opts.robustR
        if opts.niiwrite
            cd(opts.corrCVRdir);
            niftiwrite(cast(mask.*robustIR,opts.mapDatatype),'CVR_based_on_highest_r_map',opts.info.map);
            niftiwrite(cast(mask.*(robustIR.^2),opts.mapDatatype),'CVR_based_on_highest_r_r2_map',opts.info.map);
        else
            saveImageData(mask.*robustIR, opts.headers.map, opts.corrCVRdir,'CVR_based_on_highest_r_map.nii.gz', datatype);
            saveImageData(mask.*(robustIR.^2), opts.headers.map, opts.corrCVRdir,'CVR_based_on_highest_r_r2_map.nii.gz', datatype);
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
            %beta(beta > 30) = 0; beta(beta < -30) = 0;
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
            %beta(beta > 30) = 0; beta(beta < -30) = 0;
        end

        lagmatrix = lags(maxindex);
        GLM_Estimate = zeros([xx*yy*zz,1]); GLM_Estimate(coordinates,:) = beta;
        GLM_lags = zeros([xx*yy*zz,1]); GLM_lags(coordinates,:) = lagmatrix;
        GLM_lags = reshape(GLM_lags,[xx yy zz]);
        GLM_Estimate = reshape(GLM_Estimate,[xx yy zz]);

        tmp = opts.TR*(lagmatrix/opts.interp_factor);
        if min(tmp(:)) < 0; tmp = tmp + abs(min(tmp(:))); end
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

        if opts.medfilt_maps
            tmpLag = medfilt3(tmpLag);
            GLM_lags = medfilt3(GLM_lags);
        end

        switch pp
            case 1
                if opts.niiwrite
                    cd(opts.glmlagdir);
                    niftiwrite(cast(mask.*GLM_Estimate,opts.mapDatatype),'GLM_refined_probe_beta_map',opts.info.map);
                    niftiwrite(cast(mask.*tmpLag,opts.mapDatatype),'GLM_refined_probe_hemodynamic_lag_map',opts.info.map);
                    niftiwrite(cast(mask.*GLM_lags,opts.mapDatatype),'raw_GLM_refined_probe_hemodynamic_lag_map',opts.info.map);
                    niftiwrite(cast(mask.*cR2,opts.mapDatatype),'GLM_refined_probe_R2_map',opts.info.map);
                    niftiwrite(cast(mask.*cTstat,opts.mapDatatype),'GLM_refined_probe_Tstat_map',opts.info.map);
                else
                    saveImageData(mask.*GLM_Estimate, opts.headers.map, opts.glmlagdir, 'GLM_refined_probe_beta_map.nii.gz', datatype);
                    saveImageData(mask.*tmpLag, opts.headers.map, opts.glmlagdir, 'GLM_refined_probe_hemodynamic_lag_map.nii.gz', datatype);
                    saveImageData(mask.*GLM_lags, opts.headers.map,opts.glmlagdir, 'raw_GLM_refined_probe_hemodynamic_lag_map.nii.gz', datatype);
                    saveImageData(mask.*cR2, opts.headers.map, opts.glmlagdir, 'GLM_refined_probe_R2_map.nii.gz', datatype);
                    saveImageData(mask.*cTstat, opts.headers.map, opts.glmlagdir, 'GLM_refined_probe_Tstat_map.nii.gz', datatype);
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
                    niftiwrite(cast(mask.*GLM_Estimate,opts.mapDatatype),'GLM_input_probe_beta_map',opts.info.map);
                    niftiwrite(cast(mask.*tmpLag,opts.mapDatatype),'GLM_input_probe_hemodynamic_lag_map',opts.info.map);
                    niftiwrite(cast(mask.*GLM_lags,opts.mapDatatype),'raw_GLM_input_probe_hemodynamic_lag_map',opts.info.map);
                    niftiwrite(cast(mask.*cR2,opts.mapDatatype),'GLM_input_probe_R2_map',opts.info.map);
                    niftiwrite(cast(mask.*cTstat,opts.mapDatatype),'GLM_input_probe_Tstat_map',opts.info.map);
                else
                    saveImageData(mask.*GLM_Estimate, opts.headers.map, opts.glmlagdir, 'GLM_input_probe_beta_map.nii.gz', datatype);
                    saveImageData(mask.*tmpLag, opts.headers.map, opts.glmlagdir, 'GLM_input_probe_hemodynamic_lag_map.nii.gz', datatype);
                    saveImageData(mask.*GLM_lags, opts.headers.map, opts.glmlagdir, 'raw_GLM_input_probe_hemodynamic_lag_map.nii.gz', datatype);
                    saveImageData(mask.*cR2, opts.headers.map, opts.glmlagdir, 'GLM_input_probe_R2_map.nii.gz', datatype);
                    saveImageData(mask.*cTstat, opts.headers.map, opts.glmlagdir, 'GLM_input_probe_Tstat_map.nii.gz', datatype);
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
                B = wb_voxel_ts(ii,:); B(isnan(A)) = []; A(isnan(A)) = [];
                C = [ones([length(A) 1]) A'];
                regr_coef(:,ii)= C\B';
            end

            cCVR(1,coordinates) = regr_coef(end,:); %extract slope
            cCVR(cCVR > 30) = 0; cCVR(cCVR < -30) = 0; %cleanup base CVR map
            cCVR = reshape(cCVR, [xx yy zz]);
            %save lag-adjusted CVR map
            if pp == 1
                if opts.niiwrite
                    cd(opts.glmCVRdir)
                    niftiwrite(cast(mask.*cCVR,opts.mapDatatype),'GLM_refined_probe_lag_corrected_CVR_map',opts.info.map);
                else
                    saveImageData(mask.*cCVR, opts.headers.map, opts.glmCVRdir, 'GLM_refined_probe_lag_corrected_CVR_map.nii.gz', datatype);
                end
            else
                if opts.niiwrite
                    cd(opts.glmCVRdir)
                    niftiwrite(cast(mask.*cCVR,opts.mapDatatype),'GLM_input_probe_lag_corrected_CVR_map',opts.info.map);
                else
                    saveImageData(mask.*cCVR, opts.headers.map, opts.glmCVRdir, 'GLM_input_probe_lag_corrected_CVR_map.nii.gz', datatype);
                end
            end
            maps.GLM.CVR.optiReg_cCVR = cCVR;
            %calculate statistics
            SSE = zeros([1 size(wb_voxel_ts,1)]); SST = zeros([1 size(wb_voxel_ts,1)]); cT = zeros([1 size(wb_voxel_ts,1)]);
            parfor ii=1:size(wb_voxel_ts,1)
                A = shifted_regr(ii,:);
                X = wb_voxel_ts(ii,:); X(isnan(A)) = []; A(isnan(A)) = [];
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
                    niftiwrite(cast(mask.*cR2,opts.mapDatatype),'GLM_refined_probe_lag_corrected_R2_map',opts.info.map);
                    niftiwrite(cast(mask.*cTstat,opts.mapDatatype),'GLM_refined_probe_lag_corrected_Tstat_map',opts.info.map);
                else
                    saveImageData(mask.*cR2, opts.headers.map, opts.glmCVRdir, 'GLM_refined_probe_lag_corrected_R2_map.nii.gz', datatype);
                    saveImageData(mask.*cTstat, opts.headers.map, opts.glmCVRdir,'GLM_refined_probe_lag_corrected_Tstat_map.nii.gz', datatype);
                end
                maps.GLM.CVR.optiReg_cR2 = cR2;
                maps.GLM.CVR.optiReg_cSSE = cSSE;
                maps.GLM.CVR.optiReg_cTstat = cTstat;
            else
                if opts.niiwrite
                    cd(opts.glmCVRdir)
                    niftiwrite(cast(mask.*cR2,opts.mapDatatype),'GLM_input_probe_lag_corrected_R2_map',opts.info.map);
                    niftiwrite(cast(mask.*cTstat,opts.mapDatatype),'GLM_input_probe_lag_corrected_Tstat_map',opts.info.map);
                else
                    saveImageData(mask.*cR2, opts.headers.map, opts.glmCVRdir, 'GLM_input_probe_lag_corrected_R2_map.nii.gz', datatype);
                    saveImageData(mask.*cTstat, opts.headers.map, opts.glmCVRdir,'GLM_input_probe_lag_corrected_Tstat_map.nii.gz', datatype);
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
    save(fullfile(opts.resultsdir,'result_maps.mat'), 'maps');
    try
        warning('on');
        warning('on', 'MATLAB:rankDeficientMatrix');

        spmd
            warning('on', 'MATLAB:rankDeficientMatrix');
        end
    catch
    end

end

