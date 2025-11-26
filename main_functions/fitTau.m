% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% a.bhogal@umcutrecht.nl
% <fitTau: fits a convolved hemodynamic response function to extract dynamic (i.e. TAU) CVR information >
%
% (GPL header unchanged)
%
% *************************************************************************
% this function is inspited by the manuscript from Poublanc et al.:
% Measuring cerebrovascular reactivity: the dynamic response to a step hypercapnic stimulus
% DOI: 10.1038/jcbfm.2015.114
%
% data: input timeseries data (i.e. 4D BOLD MRI dataset). IMPORTANT:
% Volumes with zero values will not be weighted in the tau fit. This allows
% to weight specific aspects of the data - for example. only a rising or
% falling slope.
%
% mask: binary mask defining voxels of interest
%
% probe: usually PetCO2 or optimized regressor derived from lagCVR
% function. The probe is the 'reference' response
%
% weighting: Added a weighting such that when input data contains volumes
% (or timeseries) with zero values, those values are weighted less in the
% fit. This means that rising or falling edges can be isolated
% independently
%

function [maps, tau_fits] = fitTau(probe, data, mask, opts)
global opts
tic

ty   = class(data);
data = double(data);
mask = logical(mask);
if iscolumn(probe); probe = probe'; end

if isfield(opts,'interp_factor'); else; opts.interp_factor = 1; end
if isfield(opts,'refine_tau');    else; opts.refine_tau = 1; end
if isfield(opts,'passes');        else; opts.passes = 10; end
if isfield(opts,'win_size');      else; opts.win_size = 1; end
if isfield(opts,'max_tau');       else; opts.max_tau = 300; end
if isfield(opts,'save_unrefined');else; opts.save_unrefined = 0; end
if isfield(opts,'filloutliers');  else; opts.filloutliers = 1; end
if isfield(opts,'medfilt_maps');  else; opts.medfilt_maps = 1; end
if isfield(opts,'save_responses');else; opts.save_responses = 0; end

opts.dynamicdir = fullfile(opts.resultsdir,'tau'); mkdir(opts.dynamicdir);
savedir = opts.dynamicdir;
[xx, yy, zz, dyn] = size(data);

if opts.filloutliers
    disp('removing outliers...')
    [voxel_ts, coordinates] = grabTimeseries(data, mask);
    voxel_ts = filloutliers(voxel_ts', 'spline', 'mean');
    voxel_ts = voxel_ts';
else
    [voxel_ts, coordinates] = grabTimeseries(data, mask);
end
voxel_ts(isnan(voxel_ts)) = 0; voxel_ts(isinf(voxel_ts)) = 0;

ts = zeros([length(coordinates), opts.interp_factor*dyn]);

% ---------------------------------
% Interpolation and time vector
% ---------------------------------
if size(data,4) ~= length(probe) && opts.interp_factor > 1
    disp('interpolating data')
    t = opts.TR/opts.interp_factor : opts.TR/opts.interp_factor : size(data,4)*opts.TR;
    parfor ii = 1:length(coordinates)
        ts(ii,:) = interp(voxel_ts(ii,:), opts.interp_factor);
    end
    clear voxel_ts;
elseif size(data,4) == length(probe) && opts.interp_factor > 1
    disp('interpolating data and probe')
    t = opts.TR/opts.interp_factor : opts.TR/opts.interp_factor : size(data,4)*opts.TR;
    parfor ii = 1:length(coordinates)
        ts(ii,:) = interp(voxel_ts(ii,:), opts.interp_factor);
    end
    probe = interp(probe, opts.interp_factor);   % <-- original probe, kept in native units
    clear voxel_ts;
else
    t  = opts.TR:opts.TR:size(data,4)*opts.TR;
    ts = voxel_ts; clear voxel_ts;
end

% rescale probe ONLY for fitting the exponential model
input_probe = rescale(probe);
ts(isnan(ts)) = 0;

nr_params = 3;
b         = nan([length(coordinates), nr_params]);

probe_fft = fftinput(input_probe);
model     = @(a,t) real(ifft(ifftshift(probe_fft .* fftexponential(a(2),a(1),t)))) + a(nr_params);

% Weighted lsqnonlin settings
lb = double([-Inf, 0, -Inf]);
ub = double([ Inf, opts.max_tau, Inf]);

options = optimoptions('lsqnonlin', 'Display', 'none', ...
    'FunctionTolerance', 1.0000e-8, 'StepTolerance', 1.0000e-8, 'MaxIter', 500);

queue = createParallelProgressBar(numel(coordinates));
parfor ii = 1:length(coordinates)
    valid_idx = ~(ts(ii, :) == 0 | isnan(ts(ii, :)));
    weights   = double(valid_idx);

    weighted_residual = @(a) weights .* (ts(ii, :) - model(a, t));

    try
        b(ii, :) = lsqnonlin(weighted_residual, [range(ts(ii, :)) 20 0], lb, ub, options);
    catch
        disp(['Error voxel ', int2str(ii)]);
    end
    pause(0.00001);
    send(queue, ii);
end

b1_vec = zeros([1, xx*yy*zz]);
b2_vec = b1_vec; b3_vec = b1_vec;

b1_vec(1, coordinates) = b(:,1);
b2_vec(1, coordinates) = b(:,2);
b3_vec(1, coordinates) = b(:,3);

% refined analysis for clipped tau values
if opts.save_unrefined && opts.refine_tau
    b1_map = reshape(b1_vec, [xx yy zz]);
    b2_map = reshape(b2_vec, [xx yy zz]);
    b3_map = reshape(b3_vec, [xx yy zz]);
    savedir = opts.dynamicdir;
    saveMap(cast(mask.*b1_map,opts.mapDatatype), savedir, 'exp_scaling_unrefined',  opts.info.map, opts);
    saveMap(cast(mask.*b2_map,opts.mapDatatype), savedir, 'exp_tau_unrefined',      opts.info.map, opts);
    saveMap(cast(mask.*b3_map,opts.mapDatatype), savedir, 'exp_offset_unrefined',   opts.info.map, opts);
end

if opts.refine_tau

    b2_map = reshape(b2_vec, [xx yy zz]);
    iter   = 1;
    passes = 0;
    max_tau = max(ceil(b2_map(:)));

    while iter
        passes = passes + 1;
        max_map = zeros(size(mask));
        max_map(ceil(b2_map) == max_tau) = 1;

        % get new coordinates
        [~, newcoordinates] = grabTimeseries(data, max_map);
        clear b
        b = nan([length(newcoordinates), nr_params]);
        i = find(max_map);
        [X,Y,Z] = ind2sub(size(max_map), i); clear i
        newTS = zeros([length(X) size(data,4)*opts.interp_factor]);

        for kk=1:length(X)
            try
                X_rng = X(kk)-opts.win_size:X(kk)+opts.win_size;
                Y_rng = Y(kk)-opts.win_size:Y(kk)+opts.win_size;
                Z_rng = Z(kk)-opts.win_size:Z(kk)+opts.win_size;

                X_rng(X_rng > size(data,1)) = []; X_rng(X_rng < 1) = [];
                Y_rng(Y_rng > size(data,2)) = []; Y_rng(Y_rng < 1) = [];
                Z_rng(Z_rng > size(data,3)) = []; Z_rng(Z_rng < 1) = [];

                tmp =  data(X_rng,Y_rng,Z_rng,:);
                tmp = reshape(tmp, [size(tmp,1)*size(tmp,2)*size(tmp,3), size(tmp,4)]);
                tmp(find(tmp(:,1) == 0),:) = [];
                data(X(kk),Y(kk),Z(kk),:) = nanmean(tmp);
                newTS(kk,:) = interp(nanmean(tmp),opts.interp_factor);
            catch
                disp('error; skipping voxel')
            end
        end

        clear tmp
        parfor ii = 1:length(newcoordinates)
            valid_idx = ~(newTS(ii, :) == 0 | isnan(newTS(ii, :)));
            weights   = double(valid_idx);

            weighted_residual = @(a) weights .* (newTS(ii, :) - model(a, t));

            try
                b(ii, :) = lsqnonlin(weighted_residual, [range(newTS(ii, :)) 20 0], lb, ub, options);
            catch
                disp(['Error voxel ', int2str(ii)]);
            end
        end

        b1_vec(1, newcoordinates) = b(:,1);
        b2_vec(1, newcoordinates) = b(:,2);
        b3_vec(1, newcoordinates) = b(:,3);

        b2_map = reshape(b2_vec, [xx yy zz]);
        max_map = zeros(size(b2_map));
        max_map(ceil(b2_map) == max_tau) = 1;
        perc = 100*(nnz(max_map)/nnz(mask));
        disp([int2str(perc), ' percent of voxels have clipped tau values'])
        if perc > 2
            if passes < opts.passes
                continue;
            else
                iter = 0;
                disp('exceeded the maximum allowed passes')
                disp('... to increase passes set opts.passes to a higher value')
            end
        else
            iter = 0;
        end
    end

    b1_map = reshape(b1_vec, [xx yy zz]);
    b2_map = reshape(b2_vec, [xx yy zz]);
    b3_map = reshape(b3_vec, [xx yy zz]);
end

%% save maps and calculate stats

responseFits = zeros(length(coordinates), length(t));
b = [];
b(1,:) = b1_map(coordinates);
b(2,:) = b2_map(coordinates);
b(3,:) = b3_map(coordinates);

SSE    = zeros(1,length(coordinates));
SST    = zeros(1,length(coordinates));
r_corr = zeros(1,length(coordinates));

% ---------- NEW: tau-corrected CVR storage ----------
tau_cvr = nan(1, length(coordinates));
% FFT of ORIGINAL (possibly interpolated) probe in native units
probe_fft_orig = fftinput(probe);

parfor ii = 1:length(coordinates)
    % model fit (using rescaled probe)
    responseFits(ii,:) = model(b(:,ii), t);

    SSE(1,ii) = (norm(ts(ii,:) - responseFits(ii,:)))^2;
    SST(1,ii) = (norm(ts(ii,:)-mean(ts(ii,:))))^2;
    r_corr(1,ii) = corr(ts(ii,:)', responseFits(ii,:)');

    % ---------------------------------------------------------
    % Tau-corrected CVR:
    %   - use voxel-specific tau
    %   - convolve ORIGINAL probe with exponential dispersion
    %   - regress voxel ts against dispersed probe
    % ---------------------------------------------------------
    tau_i = b(2,ii);

    % exponential kernel with unit amplitude (tau = tau_i)
    kernel_fft = fftexponential(tau_i, 1, t);
    dispersed_probe = real(ifft(ifftshift(probe_fft_orig .* kernel_fft)));

    % same weighting strategy: ignore zeros / NaNs in ts
    valid_idx = ~(ts(ii,:) == 0 | isnan(ts(ii,:)));
    x = dispersed_probe(valid_idx)';
    y = ts(ii,valid_idx)';

    if numel(x) > 2 && any(x) && any(y)
        % linear regression with intercept: y = beta1*x + beta0
        X = [x, ones(numel(x),1)];
        beta = X \ y;
        tau_cvr(ii) = beta(1);   % slope = tau-corrected CVR
    else
        tau_cvr(ii) = NaN;
    end
end

R2  = 1 - SSE./SST;
cR2 = zeros(1, numel(mask)); 
r   = zeros(1, numel(mask));
cR2(1, coordinates) = R2; 
r(1, coordinates)   = r_corr;
cR2 = reshape(cR2, size(mask)); 
r   = reshape(r,   size(mask));

% tau-corrected CVR map (full volume)
tau_cvr_map = zeros(1, numel(mask));
tau_cvr_map(1, coordinates) = tau_cvr;
tau_cvr_map = reshape(tau_cvr_map, size(mask));

if opts.medfilt_maps
    b1_map      = medfilt3(b1_map);
    b2_map      = medfilt3(b2_map);
    b3_map      = medfilt3(b3_map);
    tau_cvr_map = medfilt3(tau_cvr_map);
end

saveMap(cast(mask.*b1_map,      opts.mapDatatype), savedir, 'signal_magnitude',     opts.info.map, opts);
saveMap(cast(mask.*b2_map,      opts.mapDatatype), savedir, 'signal_dispersion',    opts.info.map, opts);
saveMap(cast(mask.*b3_map,      opts.mapDatatype), savedir, 'signal_offset',        opts.info.map, opts);
saveMap(cast(mask.*cR2,         opts.mapDatatype), savedir, 'R2_map',               opts.info.map, opts);
saveMap(cast(mask.*r,           opts.mapDatatype), savedir, 'r_map',                opts.info.map, opts);
saveMap(cast(mask.*(r.^2),      opts.mapDatatype), savedir, 'expVariance_r2_map',   opts.info.map, opts);
saveMap(cast(mask.*tau_cvr_map, opts.mapDatatype), savedir, 'tau_corrected_CVR',    opts.info.map, opts);

if opts.save_responses
    tau_fits = zeros(numel(mask), length(t));
    tau_fits(coordinates, :) = responseFits;
    tau_fits = reshape(tau_fits, [xx, yy, zz, length(t)]);
    saveMap(cast(tau_fits, opts.tsDatatype), savedir, 'tau_fits', opts.info.ts, opts);
    maps.fitted_responses = tau_fits;
    clear tau_fits
end

maps.expHRF.signal_magnitude     = b1_map;
maps.expHRF.signal_dispersion    = b2_map;
maps.expHRF.signal_offset        = b3_map;
maps.expHRF.R2_map               = cR2;
maps.expHRF.r_map                = r;
maps.expHRF.expVariance_r2_map   = r.^2;
maps.expHRF.tau_corrected_CVR    = tau_cvr_map;

toc
disp('finished running tau analysis')
disp('...saving maps in .mat file' )

save(fullfile(opts.dynamicdir,'result_tau_maps.mat'), 'maps');
end
