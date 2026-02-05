
function [maps, BP_ref, bpData] = gaslessCVR(data, mask, refmask, nuisance, opts)
% Copyright (C) Alex A. Bhogal
% gaslessCVR: resting-state CVR via GLM of a band-passed reference regressor
% CVR index = beta for the reference regressor; nuisance (if provided) is included
% solely to remove variance, but is not of direct interest.
global opts
%% ---- Quiet common warnings (safe-guard) --------------------------------
try
    warning('off');
    warning('off', 'MATLAB:rankDeficientMatrix');
    warning('off', 'MATLAB:nearlySingularMatrix')
    warning('off', 'MATLAB:singularMatrix')
    spmd
        warning('off');
        warning('off', 'MATLAB:rankDeficientMatrix');
        warning('off', 'MATLAB:nearlySingularMatrix')
        warning('off', 'MATLAB:singularMatrix')
    end
catch
end

%% ---- Inputs & defaults --------------------------------------------------
refmask = logical(refmask);
mask    = logical(mask);

if ~isfield(opts,'fpass');         opts.fpass        = [0.000001 0.1164]; end  % Hz
if ~isfield(opts,'niiwrite');      opts.niiwrite     = 0;                  end
if ~isfield(opts,'prepNuisance');  opts.prepNuisance = 1;                  end
if ~isfield(opts,'TR'); error('opts.TR must be provided (seconds).'); end

% Nuisance orientation to T x K
if isempty(nuisance)
    np = [];
    opts.prepNuisance = 0;
else
    if size(nuisance,1) < size(nuisance,2)
        nuisance = nuisance'; % make T x K
    end
end

% Prepare nuisance (optional preprocessing)
if opts.prepNuisance
    reference  = meanTimeseries(data, refmask);      % 1 x T
    [np, ~]    = prepNuisance(nuisance, reference, opts);
    clear reference
else
    np = nuisance;
end

% Output dirs
if isempty(np) || nnz(np)==0
    opts.CVRidxdir = fullfile(opts.resultsdir,'CVRidx'); 
else
    opts.CVRidxdir = fullfile(opts.resultsdir,'CVRidx_inclNuisance');
end
if ~exist(opts.CVRidxdir,'dir'), mkdir(opts.CVRidxdir); end
if ~isfield(opts,'figdir'), opts.figdir = opts.CVRidxdir; end

%% ---- Data shapes --------------------------------------------------------
[xx,yy,zz,NT] = size(data);
data = double(data);

% Extract voxel time series
[voxels, coordinates] = grabTimeseries(data, mask);     % V x T
[ref_voxels, ~]       = grabTimeseries(data, refmask);  % Vr x T

T = size(voxels,2);   % timepoints
V = size(voxels,1);   % voxels in mask
assert(T == NT, 'Internal shape mismatch: time dimension must match.');

%% ---- Band-pass filtering of reference and voxel data --------------------
Fs   = 1/opts.TR;              % Hz
Lowf = opts.fpass(1); 
Highf= opts.fpass(2);

mean_ts = mean(ref_voxels,1);  % 1 x T

if license('test','signal_toolbox') == 1
    % Choose/compute filter
    if isfield(opts,'filter_order')
        [b,a]  = butter(opts.filter_order, 2*[Lowf, Highf]/Fs);
        BP_ref = filtfilt(b,a, mean_ts);
    else
        % Auto-select order by minimizing SSE between original and filtered
        SSres = zeros(1,4);
        BP    = zeros(4, numel(mean_ts));
        for forder = 1:4
            [bb,aa]     = butter(forder, 2*[Lowf, Highf]/Fs);
            BP(forder,:) = filtfilt(bb,aa, mean_ts);
            SSres(forder)= sum((mean_ts - BP(forder,:)).^2);
        end
        [~,I]  = min(SSres);
        [b,a]  = butter(I, 2*[Lowf, Highf]/Fs);
        BP_ref = filtfilt(b,a, mean_ts);
        opts.filter_order = I;
    end

    % Band-pass all mask voxels
    BP_V  = zeros(size(voxels));     % V x T
    parfor ii = 1:V
        BP_V(ii,:) = filtfilt(b,a, voxels(ii,:));
    end

    % Band-pass reference voxels (not strictly needed further, but kept)
    BP_rV = zeros(size(ref_voxels));
    parfor ii = 1:size(ref_voxels,1)
        BP_rV(ii,:) = filtfilt(b,a, ref_voxels(ii,:));
    end
else
    isplot = 0;
    BP_ref = bpfilt(mean_ts, Lowf, Highf, isplot); BP_ref = BP_ref';
    BP_V   = bpfilt(voxels,     Lowf, Highf, isplot);     % V x T
    BP_rV  = bpfilt(ref_voxels, Lowf, Highf, isplot);
end

%% ---- QC plot of ref regressor ------------------------------------------
figure; hold on
plot(rescale(mean_ts),'DisplayName','mean reference time-series');
plot(rescale(BP_ref)+0.2,'DisplayName','BP reference (offset)');
title('reference regressor before and after bandpass');
xlabel('image volumes'); ylabel('a.u.'); legend('show');
saveas(gcf, fullfile(opts.figdir,'reference_regressor.fig'));

%% ---- Build GLM design: intercept + (optional) nuisance + reference -----
BP_ref_scaled = rescale(BP_ref(:));       % T x 1; you can switch to zscore if desired

if isempty(np) || nnz(np)==0
    D = [ones(T,1), BP_ref_scaled];       % T x p ; p=2
    slopeIdx = 2;
else
    if size(np,1) ~= T
        error('nuisance must have T rows to match timepoints.');
    end
    D = [ones(T,1), np, BP_ref_scaled];   % T x p ; p=2+K, ref is last column
    slopeIdx = size(D,2);
end
p = size(D,2);

% Precompute XtX inverse for SE
XtX     = D.'*D;                % p x p
XtX_inv = pinv(XtX);            % robust inverse

%% ---- Solve LS for all voxels at once -----------------------------------
% We want coef as p x V
try
    coef = gather( gpuArray(D) \ gpuArray(BP_V.') );   % (T x p) \ (T x V) -> p x V
catch
    coef = D \ BP_V.';                                 % p x V
end

% Fitted values for full model
Y = D * coef;                     % (T x p)*(p x V) -> T x V

% Residuals and sums of squares
E    = BP_V.' - Y;                % T x V
SSE  = sum(E.^2, 1);              % 1 x V
muV  = mean(BP_V, 2);             % V x 1
SST  = sum( (BP_V.' - muV.').^2 , 1 );   % 1 x V
vR2  = 1 - SSE ./ max(SST, eps);  % 1 x V

% t-stat for reference beta
sigma2  = SSE ./ max(T - p, 1);                   % 1 x V
SE_beta = sqrt( XtX_inv(slopeIdx,slopeIdx) .* sigma2 ); % 1 x V
vT      = coef(slopeIdx,:) ./ max(SE_beta, eps);  % 1 x V

% CVR index = beta for reference regressor
CVRidx  = coef(slopeIdx,:);                       % 1 x V

%% ---- Pack band-passed data back into 4-D volume -------------------------
bpData = zeros(size(data));
bpData = reshape(bpData,[xx*yy*zz, NT]);
bpData(coordinates, :) = BP_V;      % V x T
bpData = reshape(bpData, size(data));

%% ---- Write maps into volumes -------------------------------------------
R2    = zeros(1, xx*yy*zz);  R2(1,coordinates)    = vR2;   R2    = reshape(R2,[xx yy zz]);
Tstat = zeros(1, xx*yy*zz);  Tstat(1,coordinates) = vT;    Tstat = reshape(Tstat,[xx yy zz]);

CVRidx_map = zeros(1, numel(mask));
CVRidx_map(1, coordinates) = CVRidx;
CVRidx_map = reshape(CVRidx_map, size(mask));

% Outlier clipping (keep your original convention)
CVRidx_map(CVRidx_map > 100)  = 0; 
CVRidx_map(CVRidx_map < -100) = 0;

% Normalize by mean ref ROI response (internally normalized CVR)
mean_ref_response = ROImean(CVRidx_map, refmask);
nCVRidx_map = CVRidx_map / mean_ref_response;

% Clip normalized outliers (same convention)
nCVRidx_map(nCVRidx_map > 100)  = 0; 
nCVRidx_map(nCVRidx_map < -100) = 0;

% Save
saveMap(nCVRidx_map, opts.CVRidxdir,'GLM_CVR_normalized_map',opts.info.map, opts);
saveMap(CVRidx_map, opts.CVRidxdir, 'GLM_CVR_non_normalized_map',opts.info.map, opts);
saveMap(R2, opts.CVRidxdir, 'GLM_CVR_R2_map',opts.info.map, opts);
saveMap(Tstat, opts.CVRidxdir, 'GLM_CVR_Tstat_map',opts.info.map, opts);

% Outputs struct
nCVRidx_map(nCVRidx_map == 0) = NaN;
CVRidx_map(CVRidx_map == 0)   = NaN;
maps.nCVRidx = nCVRidx_map;
maps.CVRidx  = CVRidx_map;
maps.R2      = R2;
maps.Tstat   = Tstat;

%% ---- Restore warnings ---------------------------------------------------
try
    warning('on');
    warning('on', 'MATLAB:rankDeficientMatrix');
    warning('on', 'MATLAB:nearlySingularMatrix')
    warning('on', 'MATLAB:singularMatrix')
    spmd
        warning('on');
        warning('on', 'MATLAB:rankDeficientMatrix');
        warning('on', 'MATLAB:nearlySingularMatrix')
        warning('on');
    end
catch
end
end
