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
% Copyright (C) Alex A. Bhogal
% sections of this code were contributed by Allen A. Champagne
% lagCVR2: calculates hemodynamic parameter maps with associated statistical maps
%
% This version fixes several correctness issues in lagCVR while preserving the
% original structure and outputs. Requested: keep the `global opts` aspect intact.

refmask = logical(refmask); 
mask    = logical(mask);
probe   = double(probe);
data    = double(data);

try
    warning('off');
    warning('off', 'MATLAB:rankDeficientMatrix');
    spmd
        warning('off', 'MATLAB:rankDeficientMatrix');
    end
catch
end

% (per user request) keep this line unchanged
global opts;

% Clean data
data(isnan(data)) = 0; 
data(isinf(data)) = 0;

% Check for stats toolbox (for PCA/licols usage flag)
if license('test','Statistics_toolbox') == 0
    opts.pca_analysis = 0;
end

maps = struct();

% Probe must be a column
if ~iscolumn(probe), probe = probe'; end

%% Nuisance handling & preparation
np = [];
norm_np = [];
idx = [];

if isempty(nuisance)
    np = [];
else
    % Ensure T x K
    if size(nuisance,1) < size(nuisance,2), nuisance = nuisance'; end
    % Interpolate & rescale (original behavior)
    np = zeros(size(nuisance,1)*opts.interp_factor, size(nuisance,2));
    for ii = 1:size(nuisance,2)
        np(:,ii) = rescale(interp(nuisance(:,ii), opts.interp_factor), -1, 1);
    end
    % Prepare nuisance with helper (original call)
    opts.motioncorr = 0.3;
    [np, ~] = prepNuisance(nuisance, probe, opts); %#ok<ASGLU>
    % Remove linearly dependent columns
    [norm_np, idx] = licols(np); %#ok<ASGLU>
end

%% Defaults
if ~isfield(opts,'gpu'),             opts.gpu = 0;                end
if ~isfield(opts,'niiwrite'),        opts.niiwrite = 0;           end
if ~isfield(opts,'plot'),            opts.plot = 0;               end
if ~isfield(opts,'prewhite'),        opts.prewhite = 0;           end
if ~isfield(opts,'interp_factor'),   opts.interp_factor = 4;      end
if ~isfield(opts,'load_probe'),      opts.load_probe = 0;         end
if ~isfield(opts,'save_rts'),        opts.save_rts = 0;           end
if ~isfield(opts,'trace_corr'),      opts.trace_corr = 1;         end
if ~isfield(opts,'refine_regressor'),opts.refine_regressor = 1;   end
if ~isfield(opts,'corr_model'),      opts.corr_model = 1;         end
if ~isfield(opts,'cvr_maps'),        opts.cvr_maps = 1;           end
if ~isfield(opts,'eff_probe'),       opts.eff_probe = 0;          end
if ~isfield(opts,'glm_model'),       opts.glm_model = 0;          end
if ~isfield(opts,'uni'),             opts.uni = 0;                end
if ~isfield(opts,'norm_regr'),       opts.norm_regr = 0;          end
if ~isfield(opts,'robust'),          opts.robust = 0;             end
if ~isfield(opts,'refine_lag'),      opts.refine_lag = 1;         end
if ~isfield(opts,'win_size'),        opts.win_size = 1;           end
if ~isfield(opts,'passes'),          opts.passes = 10;            end
if ~isfield(opts,'medfilt_maps'),    opts.medfilt_maps = 1;       end
if ~isfield(opts,'filloutliers'),    opts.filloutliers = 1;       end
if ~isfield(opts,'figdir'),          opts.figdir = opts.resultsdir; end
if opts.niiwrite
    if ~isfield(opts,'info') || ~isfield(opts.info,'rts')
        opts.info.rts = opts.info.ts;
    end
end

% Important thresholds
if ~isfield(opts,'corr_thresh'),      opts.corr_thresh = 0.7; end
if ~isfield(opts,'lowerlagthresh'),   opts.lowerlagthresh = -2; end
if ~isfield(opts,'upperlagthresh'),   opts.upperlagthresh = 2;  end

% Account for interpolation factor
opts.lowerthresh = opts.lowerlagthresh * opts.interp_factor;
opts.upperthresh = opts.upperlagthresh * opts.interp_factor;

% Correlation window (in TRs, before interp)
if ~isfield(opts,'lowlag'),  opts.lowlag  = -3; end
if ~isfield(opts,'highlag'), opts.highlag = 60; end
opts.adjlowlag  = opts.lowlag  * opts.interp_factor;
opts.adjhighlag = opts.highlag * opts.interp_factor;

% Sanity check: if everything disabled
if opts.load_probe == 0 && opts.refine_regressor == 0 && opts.trace_corr == 0
    disp('check options; stopping analysis')
    opts.glm_model = 0;
    opts.corr_model = 0;
    opts.cvr_maps = 0;
end
if opts.trace_corr == 0 && opts.robust == 1
    disp('For robust analysis, set opts.trace_corr = 1');
    disp('...continuing without creating robust maps');
    opts.robust = 0;
end

% Setup directories
if ~isfield(opts,'resultsdir'), opts.resultsdir = fullfile(pwd); end
if opts.corr_model
    opts.corrlagdir = fullfile(opts.resultsdir,'corrLAG'); if ~exist(opts.corrlagdir,'dir'), mkdir(opts.corrlagdir); end
    if opts.cvr_maps
        opts.corrCVRdir = fullfile(opts.corrlagdir,'CVR'); if ~exist(opts.corrCVRdir,'dir'), mkdir(opts.corrCVRdir); end
    end
end
if ~opts.corr_model, opts.cvr_maps = 0; end
if opts.glm_model
    opts.glmlagdir = fullfile(opts.resultsdir,'glmLAG'); if ~exist(opts.glmlagdir,'dir'), mkdir(opts.glmlagdir); end
    if opts.cvr_maps
        opts.glmCVRdir = fullfile(opts.glmlagdir,'CVR'); if ~exist(opts.glmCVRdir,'dir'), mkdir(opts.glmCVRdir); end
    end
end

% Optional image outputs
if ~isfield(opts,'robustTstat'), opts.robustTstat = 1; end
if ~isfield(opts,'robustR'),     opts.robustR = 0;   end

cd(opts.resultsdir);
datatype = 16; % NIfTI float32

[xx, yy, zz, dyn] = size(data);
orig_regr = probe;

t0 = cputime;

%% Grab coordinates / fill outliers if requested
if opts.filloutliers
    disp('removing outliers...')
    [orig_voxel_ts, coordinates] = grabTimeseries(data, mask); % V x T
    B = filloutliers(orig_voxel_ts', 'spline', 'mean');         % T x V
    B = B';
    tmp = zeros([numel(mask) size(data,4)]);
    tmp(coordinates,:) = B;
    data = reshape(tmp, size(data)); 
    clear tmp B;
else
    [orig_voxel_ts, coordinates] = grabTimeseries(data, mask); % V x T
end

% GM/reference coordinates (for potential diagnostics)
[gm_voxel_ts, gmcoordinates] = grabTimeseries(data, refmask); %#ok<NASGU>

% Prewhitening (if requested)
if opts.prewhite
    gm_voxel_ts = gm_voxel_ts';
    parfor ii = 1:size(gmcoordinates,1)
        [gm_voxel_ts(:,ii), ~, ~] = prewhiten(gm_voxel_ts(:,ii));
    end
    gm_voxel_ts = gm_voxel_ts';
    pw_voxel_ts = orig_voxel_ts';
    parfor ii = 1:size(coordinates,1)
        [pw_voxel_ts(:,ii), ~, ~] = prewhiten(pw_voxel_ts(:,ii));
    end
    pw_voxel_ts = pw_voxel_ts';
    [probe, ~, ~] = prewhiten(probe);
end

%% Optimize / load BOLD regressor (newprobe)
if opts.refine_regressor
    probename = 'final_probe.mat';
    if exist(probename,'file') && opts.load_probe
        s = load(probename); 
        newprobe = s.newprobe;
        disp('found existing probe; skipping generation of BOLD regressor');
    else
        newprobe = optimizeRegressor(probe, data, refmask, opts);
        save(fullfile(opts.resultsdir,'final_probe.mat'), 'newprobe');
        disp('Finished creating optimized regressor');
        clear s
    end
else
    newprobe = probe;
    opts.trace_corr = 0;
    % keep lags variable defined for later
    [~,lags] = xcorr(newprobe,orig_voxel_ts(1,:),'coeff'); %#ok<ASGLU>
end

%% Interpolate regressor(s)
newprobe = interp(newprobe, opts.interp_factor);   % row or column later standardized
orig_regr = interp(orig_regr, opts.interp_factor);
probe     = orig_regr;                             % keep name 'probe' for consistency

%% Build interpolated whole-brain matrix
wb_voxel_ts = zeros(length(coordinates), dyn * opts.interp_factor);
if opts.cvr_maps && opts.prewhite
    % CVR + prewhite not compatible (original note)
    opts.prewhite = 0;
    disp('prewhitening is not compatible with CVR mapping');
    disp('CVR option takes priority...');
end

if opts.prewhite
    parfor ii = 1:length(coordinates)
        wb_voxel_ts(ii,:) = interp(pw_voxel_ts(ii,:), opts.interp_factor);
    end
else
    parfor ii = 1:length(coordinates)
        wb_voxel_ts(ii,:) = interp(orig_voxel_ts(ii,:), opts.interp_factor);
    end
end

%% ------------------------------------------------------------------------
%  Correlation-based lag mapping
%% ------------------------------------------------------------------------
if opts.corr_model
    if opts.trace_corr, qq = [1 2]; else, qq = 1; end

    switch numel(qq)
        case 1
            index_map = zeros([1 numel(mask)]);
            rvec      = zeros([1, length(coordinates)]);
            index     = zeros([2 size(wb_voxel_ts,1)]);
        case 2
            index_map = zeros([2 numel(mask)]);
            rvec      = zeros([2, length(coordinates)]);
            index     = zeros([2 size(wb_voxel_ts,1)]);
    end

    for pp = qq
        switch pp
            case 1
                regr = newprobe(:)'; % optimized
                disp('performing correlation-based lag analysis using optimized regressor');
            case 2
                regr = probe(:)';    % input
                disp('performing correlation-based lag analysis using input regressor');
        end

        % Matrix cross-correlation for each voxel vs regressor
        a2 = mat2cell(wb_voxel_ts, ones(1,size(wb_voxel_ts,1)), size(wb_voxel_ts,2)); % rows -> cells
        b2 = cellfun(@(x) xcorr(x, regr, opts.adjhighlag,'coeff'), a2, 'uniformoutput', false);
        corr_probe = cell2mat(b2)';  % Lags x Voxels
        [~, lags]  = xcorr(wb_voxel_ts(1,:), regr, opts.adjhighlag,'coeff');

        % Save correlation-timeseries if requested
        if opts.save_rts
            corr_ts = zeros([xx*yy*zz size(corr_probe,1)]);
            corr_ts(coordinates,:) = corr_probe';
            corr_ts(isnan(corr_ts)) = 0;
            corr_ts = reshape(corr_ts, [xx yy zz size(corr_probe,1)]);
        end

        % Trim lag range
        idx_trim = (lags <= opts.adjlowlag) | (lags >= opts.adjhighlag);
        lags(idx_trim) = [];
        corr_probe(idx_trim,:) = [];

        % Pick max (abs) correlation per voxel
        if opts.uni
            [~, index_map(pp,coordinates)] = max(corr_probe);     % positive-only
            [~, index(pp,:)] = max(corr_probe);
        else
            [~, index_map(pp,coordinates)] = max(abs(corr_probe)); % absolute
            [~, index(pp,:)] = max(abs(corr_probe));
        end

        % r map
        parfor ii = 1:length(coordinates)
            rvec(pp,ii) = corr_probe(index(pp,ii), ii);
        end
        r_map  = zeros([xx*yy*zz 1]); 
        lag_map= zeros([1 xx*yy*zz]);
        r_map(coordinates)   = rvec(pp,:);
        r_map  = reshape(r_map, [xx yy zz]);
        lag_map(1,coordinates) = lags(1, index_map(pp,coordinates));
        lag_map = reshape(lag_map, [xx yy zz]);

        % Refinement for clipped lag values (keep original logic; fix prints)
        if opts.refine_lag
            iter   = 1;
            passes = 0;
            maxlag = max(lag_map(:));

            while iter
                passes = passes + 1;  % fixed ; to reduce console spam
                max_map = zeros(size(lag_map));
                max_map(lag_map == maxlag) = 1;

                [~, newcoordinates] = grabTimeseries(data, max_map);
                i = find(max_map);
                [X,Y,Z] = ind2sub(size(max_map), i);
                newTS = zeros([length(X) length(regr)]);

                for kk = 1:length(X)
                    try
                        X_rng = X(kk)-opts.win_size : X(kk)+opts.win_size;
                        Y_rng = Y(kk)-opts.win_size : Y(kk)+opts.win_size;
                        Z_rng = Z(kk)-opts.win_size : Z(kk)+opts.win_size;
                        % bounds
                        X_rng(X_rng > size(data,1) | X_rng < 1) = [];
                        Y_rng(Y_rng > size(data,2) | Y_rng < 1) = [];
                        Z_rng(Z_rng > size(data,3) | Z_rng < 1) = [];

                        tmp = data(X_rng, Y_rng, Z_rng, :);
                        tmp = reshape(tmp, [], size(tmp,4));
                        % keep rows with any nonzero
                        tmp = tmp(any(tmp,2), :);
                        if isempty(tmp), continue; end
                        data(X(kk),Y(kk),Z(kk),:) = mean(tmp, 1);
                        newTS(kk,:) = interp(mean(tmp,1), opts.interp_factor);
                    catch
                        % skip voxel
                    end
                end

                % Recompute correlations for refined TS
                a2 = mat2cell(newTS, ones(1,size(newTS,1)), size(newTS,2));
                b2 = cellfun(@(x) xcorr(x, regr, opts.adjhighlag, 'coeff'), a2, 'uniformoutput', false);
                corr_probe = cell2mat(b2)'; 
                corr_probe(idx_trim,:) = []; % same trim

                % Find new best lags
                rvec2  = zeros([1, length(X)]);
                index2 = zeros([1, size(newTS,1)]);
                if opts.uni
                    [~, index2] = max(corr_probe);
                else
                    [~, index2] = max(abs(corr_probe));
                end
                parfor ii = 1:length(X)
                    rvec2(1,ii) = corr_probe(index2(1,ii), ii);
                end

                r_map = r_map(:);
                lag_map = lag_map(:);
                r_map(newcoordinates)   = rvec2;
                r_map = reshape(r_map, [xx yy zz]);
                lag_map(newcoordinates) = lags(1, index2);
                lag_map = reshape(lag_map, [xx yy zz]);
                index_map(pp, newcoordinates) = index2;

                perc = 100 * numel(index2) / numel(coordinates);
                disp([int2str(perc), ' percent of voxels have clipped lag values']);
                if perc > 2
                    if passes < opts.passes
                        continue;
                    else
                        iter = 0;
                        disp('exceeded the maximum allowed passes');
                        disp('... to increase passes set opts.passes to a higher value');
                    end
                else
                    iter = 0;
                end
            end
        end

        % Save correlation-based lag & r maps
        tmpLag = opts.TR * (lag_map / opts.interp_factor) .* mask;
        if min(tmpLag(:)) < 0, tmpLag = tmpLag + abs(min(tmpLag(:))); end
        if opts.medfilt_maps
            tmpLag  = medfilt3(tmpLag);
            lag_map = medfilt3(lag_map);
        end

        savedir = opts.corrlagdir;
        switch pp
            case 1
                saveMap(mask.*tmpLag, savedir, 'hemodynamic_lag_map_refined_probe', opts.info.map, opts);
                saveMap(mask.*lag_map, savedir, 'raw_hemodynamic_lag_map_refined_probe', opts.info.map, opts);
                saveMap(mask.*r_map, savedir, 'r_map_refined_probe', opts.info.map, opts);

                maps.XCORR.lag_opti       = tmpLag;
                maps.XCORR.uncorrlag_opti = lag_map;
                maps.XCORR.r_opti         = r_map;

                if opts.save_rts
                    if opts.niiwrite
                        cd(opts.corrlagdir);
                        disp('saving correlation timeseries based on optimized regressor');
                        niftiwrite(cast(corr_ts,opts.mapDatatype),'correlation_timeseries_refined_probe',opts.info.rts);
                    else
                        disp('saving correlation timeseries based on optimized regressor');
                        saveImageData(corr_ts, opts.headers.ts, opts.corrlagdir,  'correlation_timeseries_refined_probe.nii.gz', datatype);
                    end
                end
                if opts.robust
                    LAG(1,:,:,:) = mask.*tmpLag; %#ok<AGROW>
                    RL(1,:,:,:)  = mask.*r_map;  %#ok<AGROW>
                end

            case 2
                saveMap(mask.*tmpLag, savedir, 'hemodynamic_lag_map_input_probe', opts.info.map, opts);
                saveMap(mask.*lag_map, savedir, 'raw_hemodynamic_lag_map_input_probe', opts.info.map, opts);
                saveMap(mask.*r_map, savedir, 'r_map_input_probe', opts.info.map, opts);

                maps.XCORR.lag_input       = tmpLag;
                maps.XCORR.uncorrlag_input = lag_map;
                maps.XCORR.r_input         = r_map;

                if opts.save_rts
                    if opts.niiwrite
                        cd(opts.corrlagdir);
                        niftiwrite(cast(corr_ts,opts.mapDatatype),'correlation_timeseries',opts.info.rts);
                    else
                        saveImageData(corr_ts, opts.headers.ts, opts.corrlagdir, 'correlation_timeseries.nii.gz', datatype);
                    end
                end
                if opts.robust
                    LAG(2,:,:,:) = mask.*tmpLag; %#ok<AGROW>
                    RL(2,:,:,:)  = mask.*r_map;  %#ok<AGROW>
                end
        end
        clear corr_ts r_map lag_map tmpLag a2 b2 corr_probe
    end

    index = index_map(:, coordinates); %#ok<NASGU>
    disp(['Lag, regressor and r_maps were created in: ',int2str((cputime-t0)/60),' minutes'])
end

%% Robust lag (if both input & optimized were correlated)
if opts.trace_corr && opts.robust
    disp('Calculating robust lag map (r-weighted)');
    [~,IR] = max(abs(RL),[],1,'omitnan'); IR = squeeze(IR);
    LAG2 = reshape(LAG, 2, xx*yy*zz);
    IR   = reshape(IR, 1, []);
    robustIR = zeros(size(IR));
    for ii = 1:length(IR)
        robustIR(1,ii) = LAG2(IR(1,ii), ii);
    end
    robustIR = reshape(robustIR, size(mask));
    saveMap(logical(mask).*robustIR, opts.corrlagdir, 'robust_hemodymic_lag_map_r', opts.info.map, opts);
    maps.XCORR.robust_lag_r = robustIR;
    clear LAG RL LAG2 IR robustIR
else
    clear LAG RL
end

%% ------------------------------------------------------------------------
%  CVR maps (base and lag-corrected) from correlation-based lags
%% ------------------------------------------------------------------------
if opts.cvr_maps
    if opts.trace_corr, qq = [1 2]; else, qq = 1; end
    for pp = qq
        switch pp
            case 1
                disp('Generating base maps using probe trace');
                regr = orig_regr(:)'; % row
            case 2
                % Effective probe = linear fit from BOLD regressor to probe
                rs_newprobe = rescale(newprobe, 0.0001, 1);
                coef = glmfit(rs_newprobe, probe);
                eff_probe = coef(2,1)*rs_newprobe + coef(1,1);
                if ~isempty(opts.figdir)
                    f = figure('visible','off'); 
                    subplot(3,1,1); plot(probe,'b'); title('initial probe'); ylabel('mmHg')
                    subplot(3,1,2); plot(newprobe);  title('BOLD regressor'); ylabel('% BOLD')
                    subplot(3,1,3); plot(probe,'b'); hold on; plot(eff_probe,'r'); title('effective probe'); ylabel('mmHg')
                    saveas(f, fullfile(opts.figdir,'effective_probes.fig')); close(f);
                end
                disp('Generating base maps using effective probe trace');
                regr = eff_probe(:)'; % row
        end

        % --- Base CVR (no lag correction) ---
        A = regr;
        C = [ones([length(A) 1]) A(:)];  % T x 2
        regr_coef = C \ wb_voxel_ts';    % 2 x V

        bCVR = zeros([1 xx*yy*zz]); 
        bR2  = zeros([1 xx*yy*zz]); 
        bSSE = zeros([1 xx*yy*zz]); 
        bTstat=zeros([1 xx*yy*zz]);

        bCVR(1,coordinates) = regr_coef(2,:); % slope
        bCVR(bCVR > 30)  = 0; 
        bCVR(bCVR < -30) = 0;
        bCVR = reshape(bCVR, [xx yy zz]);

        % fitted values
        X = wb_voxel_ts;     % V x T
        Y = A(:) * regr_coef(2,:) + ones(length(A),1) * regr_coef(1,:);
        
        SSE = zeros([1 size(X,1)]); 
        SST = zeros([1 size(X,1)]); 
        STDEVr = zeros([1 size(X,1)]);
        parfor ii = 1:size(X,1)
            STDEVr(1,ii) = nanstd(X(ii,:) - Y(:,ii)'); % as in original code
            SSE(1,ii)    = (norm(X(ii,:) - Y(:,ii)'))^2;
            SST(1,ii)    = (norm(X(ii,:) - mean(X(ii,:))))^2;
        end
        bT = regr_coef(2,:) ./ STDEVr(1,:);
        R2 = 1 - SSE ./ SST;

        bR2(1, coordinates)   = R2;  bR2  = reshape(bR2,  [xx yy zz]);
        bSSE(1, coordinates)  = SSE; bSSE = reshape(bSSE, [xx yy zz]);
        bTstat(1,coordinates) = bT;  bTstat=reshape(bTstat,[xx yy zz]);

        % Save (base)
        if pp == 1
            saveMap(mask.*bCVR, opts.corrCVRdir, 'lin_regr_CVR_map',       opts.info.map, opts);
            saveMap(mask.*bR2, opts.corrCVRdir, 'lin_regr_CVR_R2_map',    opts.info.map, opts);
            saveMap(mask.*bTstat, opts.corrCVRdir, 'lin_regr_CVR_Tstat_map', opts.info.map, opts);
            maps.XCORR.CVR.bCVR   = bCVR;  maps.XCORR.CVR.bR2   = bR2;
            maps.XCORR.CVR.bSSE   = bSSE;  maps.XCORR.CVR.bTstat= bTstat;
            CVR(1,:,:,:)   = mask.*bCVR;   
            TSTAT(1,:,:,:) = mask.*bTstat; 
            RC(1,:,:,:)    = mask.*bR2;    
        else
            saveMap(mask.*bCVR, opts.corrCVRdir, 'lin_regr_CVR_effective_probe_map',       opts.info.map, opts);
            saveMap(mask.*bR2, opts.corrCVRdir, 'lin_regr_CVR_effective_probe_R2_map',    opts.info.map, opts);
            saveMap(mask.*bTstat, opts.corrCVRdir, 'lin_regr_CVR_effective_probe_Tstat_map', opts.info.map, opts);
            maps.XCORR.CVR.bCVR_eff   = bCVR;  maps.XCORR.CVR.bR2_eff   = bR2;
            maps.XCORR.CVR.bSSE_eff   = bSSE;  maps.XCORR.CVR.bTstat_eff= bTstat;
            CVR(2,:,:,:)   = mask.*bCVR;   
            TSTAT(2,:,:,:) = mask.*bTstat; 
            RC(2,:,:,:)    = mask.*bR2;    
        end

        % --- Lag-corrected CVR (shift per voxel using correlation-derived indices) ---
        disp('Generating lag-adjusted maps');
        cCVR  = zeros([1 xx*yy*zz]); 
        cR2   = zeros([1 xx*yy*zz]); 
        cSSE  = zeros([1 xx*yy*zz]); 
        cTstat= zeros([1 xx*yy*zz]);

        % Use indices from correlation stage (index_map)
        idx_use = index_map( min(pp,size(index_map,1)) , coordinates ); % pick row 1 for optimized, 2 for input if present

        shifted_regr = NaN([length(coordinates), dyn * opts.interp_factor]);
        regr_row = regr(:)'; % row
        for ii = 1:length(coordinates)
            s = idx_use(ii);
            corr_regr = circshift(regr_row, s);
            if s > 0
                corr_regr(1:s) = regr_row(1);
            elseif s < 0
                corr_regr(end+s+1:end) = regr_row(end);
            end
            shifted_regr(ii,:) = corr_regr;
        end

        regr_coef = zeros([2 length(coordinates)]);
        parfor ii = 1:length(coordinates)
            A = shifted_regr(ii,:);
            B = wb_voxel_ts(ii,:);
            keep = ~isnan(A);
            A = A(keep);
            B = B(keep);
            C = [ones([length(A) 1]) A(:)];
            regr_coef(:,ii) = C \ B(:);
        end

        cCVR(1,coordinates) = regr_coef(2,:); % slope
        cCVR(cCVR > 30)  = 0; 
        cCVR(cCVR < -30) = 0; 
        cCVR = reshape(cCVR, [xx yy zz]);

        % Stats for lag-corrected
        SSE = zeros([1 size(wb_voxel_ts,1)]); 
        SST = zeros([1 size(wb_voxel_ts,1)]); 
        cT  = zeros([1 size(wb_voxel_ts,1)]);
        parfor ii = 1:size(wb_voxel_ts,1)
            A = shifted_regr(ii,:);
            Xv = wb_voxel_ts(ii,:);
            keep = ~isnan(A);
            A = A(keep);
            Xv= Xv(keep);
            Y = regr_coef(2,ii)*A(:) + regr_coef(1,ii);
            cT(1,ii) = regr_coef(2,ii) ./ nanstd(Xv - Y(:)'); % original style t~
            SSE(1,ii)= (norm(Xv - Y(:)'))^2;
            SST(1,ii)= (norm(Xv - mean(Xv)))^2;
        end
        R2 = 1 - SSE ./ SST;

        cR2(1,coordinates)   = R2;  cR2  = reshape(cR2, [xx yy zz]);
        cSSE(1,coordinates)  = SSE; cSSE = reshape(cSSE,[xx yy zz]);
        cTstat(1,coordinates)= cT;  cTstat=reshape(cTstat,[xx yy zz]);

        if pp == 1
            saveMap(mask.*cCVR, opts.corrCVRdir, 'lag_corrected_CVR_map',       opts.info.map, opts);
            saveMap(mask.*cR2, opts.corrCVRdir, 'lag_corrected_CVR_R2_map',    opts.info.map, opts);
            saveMap(mask.*cTstat, opts.corrCVRdir, 'lag_corrected_CVR_Tstat_map', opts.info.map, opts);
            maps.XCORR.CVR.cCVR   = cCVR;  maps.XCORR.CVR.cR2   = cR2;
            maps.XCORR.CVR.cSSE   = cSSE;  maps.XCORR.CVR.cTstat= cTstat;

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
                CVR(3,:,:,:)   = mask.*cCVR;   
                TSTAT(3,:,:,:) = mask.*cTstat; 
                RC(3,:,:,:)    = mask.*cR2;    
            end
        else
            saveMap(mask.*cCVR, opts.corrCVRdir, 'lag_corrected_effective_CVR_map',       opts.info.map, opts);
            saveMap(mask.*cR2, opts.corrCVRdir, 'lag_corrected_effective_CVR_R2_map',    opts.info.map, opts);
            saveMap(mask.*cTstat,opts.corrCVRdir, 'lag_corrected_effective_CVR_Tstat_map', opts.info.map, opts);
            maps.XCORR.CVR.cCVR_eff   = cCVR;  maps.XCORR.CVR.cR2_eff   = cR2;
            maps.XCORR.CVR.cSSE_eff   = cSSE;  maps.XCORR.CVR.cTstat_eff= cTstat;

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
                CVR(4,:,:,:)   = mask.*cCVR;   
                TSTAT(4,:,:,:) = mask.*cTstat; 
                RC(4,:,:,:)    = mask.*cR2;    
            end
        end
        disp('Stimulus response data is saved');
    end
end

%% Robust CVR maps (if both traces + robust requested)
if opts.cvr_maps && opts.eff_probe && opts.trace_corr && opts.robust
    disp('Calculating robust response maps');
    [~,IT] = max(abs(TSTAT),[],1,'omitnan'); IT = squeeze(IT);
    [~,IR] = max(abs(RC),[],1,'omitnan');    IR = squeeze(IR);
    CVR2 = reshape(CVR, size(CVR,1), xx*yy*zz);
    IT = reshape(IT,1,[]);
    IR = reshape(IR,1,[]);
    robustIT = zeros(size(IT));
    robustIR = zeros(size(IR));
    for ii = 1:length(IT)
        robustIT(1,ii) = CVR2(IT(1,ii), ii);
        robustIR(1,ii) = CVR2(IR(1,ii), ii);
    end
    robustIT = reshape(robustIT, size(mask));
    robustIR = reshape(robustIR, size(mask));

    if opts.robustTstat
        if opts.niiwrite
            cd(opts.corrCVRdir);
            niftiwrite(cast(mask.*robustIT, opts.mapDatatype),'CVR_based_on_highest_Tstat_map',opts.info.map);
        else
            saveImageData(mask.*robustIT, opts.headers.map, opts.corrCVRdir,'CVR_based_on_highest_Tstat_map.nii.gz', datatype);
        end
        maps.XCORR.CVR.robustCVR_TSTAT = robustIT;
    end
    if opts.robustR
        saveMap(mask.*robustIR, opts.corrCVRdir, 'CVR_based_on_highest_r_map',     opts.info.map, opts);
        saveMap(mask.*(robustIR.^2), opts.corrCVRdir, 'CVR_based_on_highest_r_r2_map',  opts.info.map, opts);
        maps.XCORR.CVR.robustCVR_R = robustIR;
    end
    clear TSTAT RC CVR CVR2 robustIT robustIR
else
    clear TSTAT RC CVR
end

%% ------------------------------------------------------------------------
%  GLM-based lag mapping (optional)
%% ------------------------------------------------------------------------
if opts.glm_model
    % Lags from autocorrelation of newprobe to set bounds (as original pattern)
    [~,lags] = xcorr(newprobe, newprobe, opts.adjhighlag,'coeff');
    idx_trim = (lags <= opts.adjlowlag) | (lags >= opts.adjhighlag);
    lags(idx_trim) = [];

    q0 = cputime;
    if opts.trace_corr, qq = [1 2]; else, qq = 1; end
    for pp = qq
        switch pp
            case 1
                disp('performing GLM-based lag analysis using OPTIMIZED regressor');
                regr = newprobe(:)'; % row
            case 2
                disp('performing GLM-based lag analysis using INPUT regressor');
                regr = probe(:)';    % row
        end

        % Build lagged regressor matrix
        Tfull = dyn * opts.interp_factor;
        regr_matrix = zeros(numel(lags), Tfull);
        for ii = 1:numel(lags)
            s = lags(ii);
            rr = circshift(regr, s);
            if s > 0
                rr(1:s) = regr(1);
            elseif s < 0
                rr(end+s+1:end) = regr(end);
            end
            regr_matrix(ii,:) = rr;
        end

        % Regress at each lag
        if isempty(norm_np) || nnz(norm_np) == 0
            regr_coef = zeros([numel(lags), 2, length(coordinates)]);
            if opts.gpu
                wbG = gpuArray(wb_voxel_ts);
                for ii = 1:size(regr_matrix,1)
                    A = regr_matrix(ii,:);
                    keep = ~isnan(A);
                    C = [ones([nnz(keep) 1]) A(keep)'];
                    regr_coef(ii,:,:) = gather(C \ wbG(:,keep)');
                end
            else
                parfor ii = 1:size(regr_matrix,1)
                    A = regr_matrix(ii,:);
                    keep = ~isnan(A);
                    C = [ones([nnz(keep) 1]) A(keep)'];
                    regr_coef(ii,:,:) = C \ wb_voxel_ts(:,keep)';
                end
            end

            % Select best lag per voxel by R2
            maxindex = zeros([1 size(wb_voxel_ts,1)]);
            beta     = zeros([1 size(wb_voxel_ts,1)]);
            rsquared = zeros([1 size(wb_voxel_ts,1)]);

            for ii = 1:size(wb_voxel_ts,1)
                tcoef = squeeze(regr_coef(:,:,ii)); % L x 2
                Xv    = wb_voxel_ts(ii,:);         % 1 x T
                R2 = zeros(1,size(regr_matrix,1));
                for jj = 1:size(regr_matrix,1)
                    keep = ~isnan(regr_matrix(jj,:));
                    tY   = tcoef(jj,2).*regr_matrix(jj,keep) + tcoef(jj,1);
                    R2(jj) = corr(tY(:), Xv(keep)')^2;
                end
                [M,I] = max(R2);
                maxindex(1,ii) = I;
                rsquared(1,ii) = M;
                beta(1,ii)     = tcoef(I,2); % slope is col 2 (no nuisance)
            end
        else
            % With nuisance: design = [1 nuis... reg]
            regr_coef = zeros([numel(lags), (size(norm_np,2)+2), length(coordinates)]);
            if opts.gpu
                wbG = gpuArray(wb_voxel_ts);
                for ii = 1:size(regr_matrix,1)
                    A = regr_matrix(ii,:);
                    keep = ~isnan(A);
                    C = [ones([nnz(keep) 1]) norm_np(keep,:) A(keep)'];
                    regr_coef(ii,:,:) = gather(C \ wbG(:,keep)');
                end
            else
                parfor ii = 1:size(regr_matrix,1)
                    A = regr_matrix(ii,:);
                    keep = ~isnan(A);
                    C = [ones([nnz(keep) 1]) norm_np(keep,:) A(keep)'];
                    regr_coef(ii,:,:) = C \ wb_voxel_ts(:,keep)';
                end
            end

            % Select best lag per voxel by R2 with full fitted Y
            maxindex = zeros([1 size(wb_voxel_ts,1)]);
            beta     = zeros([1 size(wb_voxel_ts,1)]);
            rsquared = zeros([1 size(wb_voxel_ts,1)]);
            for ii = 1:size(wb_voxel_ts,1)
                tcoef = squeeze(regr_coef(:,:,ii)); % L x (2+K)
                Xv    = wb_voxel_ts(ii,:);
                R2 = zeros(1,size(regr_matrix,1));
                for jj = 1:size(regr_matrix,1)
                    keep = ~isnan(regr_matrix(jj,:));
                    Xdes = [ones(nnz(keep),1), norm_np(keep,:), regr_matrix(jj,keep)'];
                    tY   = Xdes * tcoef(jj,:)';
                    R2(jj) = corr(tY(:), Xv(keep)')^2;
                end
                [M,I] = max(R2);
                maxindex(1,ii) = I;
                rsquared(1,ii) = M;
                % regressor β is LAST column when nuisance present
                beta(1,ii)     = tcoef(I, end);
            end
        end

        % Compose outputs
        lagmatrix   = lags(maxindex);
        GLM_Estimate= zeros([xx*yy*zz,1]); GLM_Estimate(coordinates,:) = beta;
        GLM_lags    = zeros([xx*yy*zz,1]); GLM_lags(coordinates,:) = lagmatrix;
        GLM_lags    = reshape(GLM_lags,[xx yy zz]);
        GLM_Estimate= reshape(GLM_Estimate,[xx yy zz]);

        tmp = opts.TR * (lagmatrix / opts.interp_factor);
        if min(tmp(:)) < 0, tmp = tmp + abs(min(tmp(:))); end
        tmpLag = zeros([xx*yy*zz,1]);  tmpLag(coordinates,:) = tmp;
        tmpLag = reshape(tmpLag,[xx yy zz]);

        % Stats at selected lag per voxel
        SSE = zeros([1 size(wb_voxel_ts,1)]); 
        SST = zeros([1 size(wb_voxel_ts,1)]);
        cT  = zeros([1 size(wb_voxel_ts,1)]);
        for ii = 1:size(wb_voxel_ts,1)
            jj = maxindex(ii);
            A = regr_matrix(jj,:);
            keep = ~isnan(A);
            Xv = wb_voxel_ts(ii,keep);
            if isempty(norm_np) || nnz(norm_np) == 0
                tcoef = squeeze(regr_coef(jj,:,ii)); % [intercept, beta]
                Y = tcoef(2)*A(keep) + tcoef(1);
                cT(1,ii) = tcoef(2) ./ nanstd(Xv - Y);
            else
                tcoef = squeeze(regr_coef(jj,:,ii));
                Xdes  = [ones(nnz(keep),1), norm_np(keep,:), A(keep)'];
                Y     = Xdes * tcoef(:);
                cT(1,ii) = tcoef(end) ./ nanstd(Xv - Y(:)');
            end
            SSE(1,ii) = (norm(Xv - Y(:)'))^2;
            SST(1,ii) = (norm(Xv - mean(Xv)))^2;
        end
        R2 = 1 - SSE ./ SST;

        cR2 = zeros([1 xx*yy*zz]); cSSE = zeros([1 xx*yy*zz]); cTstat = zeros([1 xx*yy*zz]);
        cR2(1, coordinates)   = R2;   cR2   = reshape(cR2, [xx yy zz]);
        cSSE(1, coordinates)  = SSE;  cSSE  = reshape(cSSE,[xx yy zz]);
        cTstat(1,coordinates) = cT;   cTstat= reshape(cTstat,[xx yy zz]);

        if opts.medfilt_maps
            tmpLag   = medfilt3(tmpLag);
            GLM_lags = medfilt3(GLM_lags);
        end

        % Save GLM outputs
        if pp == 1
            saveMap(mask.*GLM_Estimate, opts.glmlagdir, 'GLM_refined_probe_beta_map',              opts.info.map, opts);
            saveMap(mask.*tmpLag, opts.glmlagdir, 'GLM_refined_probe_hemodynamic_lag_map',  opts.info.map, opts);
            saveMap(mask.*GLM_lags, opts.glmlagdir, 'raw_GLM_refined_probe_hemodynamic_lag_map', opts.info.map, opts);
            saveMap(mask.*cR2, opts.glmlagdir, 'GLM_refined_probe_R2_map',                opts.info.map, opts);
            saveMap(mask.*cTstat, opts.glmlagdir, 'GLM_refined_probe_Tstat_map',             opts.info.map, opts);

            maps.GLM.optiReg_ES            = GLM_Estimate;
            maps.GLM.optiReg_lags          = tmpLag;
            maps.GLM.uncor_optiReg_lags    = GLM_lags;
            maps.GLM.optiReg_R2            = cR2;
            maps.GLM.optiReg_SSE           = cSSE;
            maps.GLM.optiReg_Tstat         = cTstat;
        else
            saveMap(mask.*GLM_Estimate, opts.glmlagdir, 'GLM_input_probe_beta_map',               opts.info.map, opts);
            saveMap(mask.*tmpLag, opts.glmlagdir, 'GLM_input_probe_hemodynamic_lag_map',     opts.info.map, opts);
            saveMap(mask.*GLM_lags, opts.glmlagdir, 'raw_GLM_input_probe_hemodynamic_lag_map', opts.info.map, opts);
            saveMap(mask.*cR2, opts.glmlagdir, 'GLM_input_probe_R2_map',                   opts.info.map, opts);
            saveMap(mask.*cTstat, opts.glmlagdir, 'GLM_input_probe_Tstat_map',                opts.info.map, opts);

            maps.GLM.inputReg_ES           = GLM_Estimate;
            maps.GLM.inputReg_lags         = tmpLag;
            maps.GLM.uncor_inputReg_lags   = GLM_lags;
            maps.GLM.inputReg_R2           = cR2;
            maps.GLM.inputReg_SSE          = cSSE;
            maps.GLM.inputReg_Tstat        = cTstat;
        end

        % Optional: GLM-based lag-corrected CVR maps
        if opts.cvr_maps
            disp('Generating lag-adjusted CVR maps based on GLM analysis');
            cCVR = zeros([1 xx*yy*zz]);
            shifted_regr = NaN([length(coordinates), dyn*opts.interp_factor]);

            % Use probe here (as in original)
            regr_row = probe(:)'; 
            for ii = 1:length(coordinates)
                s = maxindex(ii); % voxel-specific lag choice
                rr = circshift(regr_row, s);
                if s > 0
                    rr(1:s) = regr_row(1);
                elseif s < 0
                    rr(end+s+1:end) = regr_row(end);
                end
                shifted_regr(ii,:) = rr;
            end

            regr_coef = zeros([2 length(coordinates)]);
            parfor ii = 1:length(coordinates)
                A = shifted_regr(ii,:);
                B = wb_voxel_ts(ii,:);
                keep = ~isnan(A);
                A = A(keep); B = B(keep);
                C = [ones([length(A) 1]) A(:)];
                regr_coef(:,ii) = C \ B(:);
            end

            cCVR(1,coordinates) = regr_coef(2,:);
            cCVR(cCVR > 30)  = 0; 
            cCVR(cCVR < -30) = 0;
            cCVR = reshape(cCVR, [xx yy zz]);

            if pp == 1
                saveMap(mask.*cCVR, opts.glmCVRdir, 'GLM_refined_probe_lag_corrected_CVR_map', opts.info.map, opts);
            else
                saveMap(mask.*cCVR, opts.glmCVRdir, 'GLM_input_probe_lag_corrected_CVR_map',   opts.info.map, opts);
            end
            maps.GLM.CVR.optiReg_cCVR = cCVR;

            % Stats
            SSE = zeros([1 size(wb_voxel_ts,1)]); 
            SST = zeros([1 size(wb_voxel_ts,1)]); 
            cT  = zeros([1 size(wb_voxel_ts,1)]);
            parfor ii = 1:size(wb_voxel_ts,1)
                A = shifted_regr(ii,:);
                Xv = wb_voxel_ts(ii,:);
                keep = ~isnan(A);
                A = A(keep); Xv = Xv(keep);
                Y = regr_coef(2,ii)*A(:) + regr_coef(1,ii);
                cT(1,ii) = regr_coef(2,ii) ./ nanstd(Xv - Y(:)');
                SSE(1,ii)= (norm(Xv - Y(:)'))^2;
                SST(1,ii)= (norm(Xv - mean(Xv)))^2;
            end
            R2 = 1 - SSE ./ SST;

            cR2 = zeros([1 xx*yy*zz]); cSSE = zeros([1 xx*yy*zz]);  cTstat = zeros([1 xx*yy*zz]);
            cR2(1, coordinates)   = R2;   cR2   = reshape(cR2, [xx yy zz]);
            cSSE(1, coordinates)  = SSE;  cSSE  = reshape(cSSE,[xx yy zz]);
            cTstat(1,coordinates) = cT;   cTstat= reshape(cTstat,[xx yy zz]);

            if pp == 1
                saveMap(mask.*cR2, opts.glmCVRdir, 'GLM_refined_probe_lag_corrected_R2_map',   opts.info.map, opts);
                saveMap(mask.*cTstat, opts.glmCVRdir, 'GLM_refined_probe_lag_corrected_Tstat_map',opts.info.map, opts);
                maps.GLM.CVR.optiReg_cR2   = cR2;
                maps.GLM.CVR.optiReg_cSSE  = cSSE;
                maps.GLM.CVR.optiReg_cTstat= cTstat;
            else
                saveMap(mask.*cR2, opts.glmCVRdir, 'GLM_input_probe_lag_corrected_R2_map',     opts.info.map, opts);
                saveMap(mask.*cTstat, opts.glmCVRdir, 'GLM_input_probe_lag_corrected_Tstat_map',  opts.info.map, opts);
                maps.GLM.CVR.inputReg_cR2   = cR2;
                maps.GLM.CVR.inputReg_cSSE  = cSSE;
                maps.GLM.CVR.inputReg_cTstat= cTstat;
            end
        end
    end
    disp(['finished running GLM analysis in: ',int2str((cputime-q0)/60),' minutes'])
    disp('saving maps in .mat file');
    maps.newprobe = newprobe;
    save(fullfile(opts.resultsdir,'result_lagCVR_maps.mat'), 'maps');
    try
        warning('on');
        warning('on', 'MATLAB:rankDeficientMatrix');
        spmd
            warning('on', 'MATLAB:rankDeficientMatrix');
        end
    catch
    end
end

end
