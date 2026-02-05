% =====================================================================
% lagVB.m  –  Bayesian (grid) lag mapping for CVR
% Copyright (C) Alex A. Bhogal, 2025, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
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
% ---------------------------------------------------------------------
% Discrete‑lag Bayesian approach based on Hayes et al. (bioRxiv 2024‑02‑06).
% ---------------------------------------------------------------------
%   • Shifts a CO₂ regressor across a lag grid
%   • Adds polynomial drift + (optional) motion nuisance regressors.
%   • Computes evidence for each lag with a closed‑form Gaussian model.
%   • Picks the lag that maximises evidence per voxel.
%   • Computes softmax-weighted posterior mean and variance of lag.
% ---------------------------------------------------------------------
% Usage
%   maps = lagVB(mask, data, probe, nuis, opts)
% ---------------------------------------------------------------------
% Outputs (single‑precision maps)
%   maps.delay    – best lag [s]
%   maps.beta     – β at best lag
%   maps.R2       – R² of full model at best lag
%   maps.sigma2   – residual variance
%   maps.postLag  – posterior mean lag [s]
%   maps.postVar  – posterior variance of lag [s²]
%% ---------- defaults ---------------------------------------------------
% Prior configuration
% - opts.priorVar controls the prior variance on β (CVR amplitude). A large value weakly constrains β, useful in GM.
% - opts.delayPriorSD sets the standard deviation (in seconds) of the Gaussian delay prior. Smaller values bias toward near-zero lags.
% - For tissue-adaptive priors, you can provide:
%     opts.priorVar_GM: variance for gray matter
%     opts.priorVar_WM: variance for white matter
%     opts.delayPriorSD_GM: delay SD for gray matter
%     opts.delayPriorSD_WM: delay SD for white matter

%priorVar = 1000 → very weak constraint, wide prior, good for flexible GM fits.
%priorVar = 100 → moderate constraint, good balance between prior and data.
%priorVar = 10 → tight constraint, useful in low SNR voxels or WM to dampen β inflation.

% ---------------------------------------------------------------------
% 🧠 Prior Variance (`priorVar`) Summary
%
% Units:
%   priorVar has units of (output units)^2 / (input units)^2
%
%   Example:
%     If BOLD is in % signal change, and probe is in mmHg,
%     then priorVar is in (%BOLD/mmHg)^2
%
% Guidelines:
%
%   Scenario                        | Units             | Example β | Recommended priorVar
%  ----------------------------------------------------------------------
%   CO₂ in mmHg, BOLD in %          | (%BOLD/mmHg)^2    | 0.5       | 0.25 to 2.0
%   GM signal probe (z-scored)      | unitless          | ~0.2–0.6  | 0.1 to 0.5
%   Weak signal (resting-state)     | same              | small     | 0.01 to 0.1
%   Strong stimulus (e.g. breath-hold)| same            | large     | 1.0 to 4.0
%
% Notes:
%   • Smaller priorVar = tighter prior = more regularization
%   • Larger priorVar = looser prior = lets data dominate
%   • Use tighter prior in WM or with noisy regressors
% ---------------------------------------------------------------------
% ⏱️ Delay Prior Standard Deviation (`opts.delayPriorSD`)
%
% Units:
%   The delayPriorSD parameter is specified in **seconds**.
%
% Description:
%   Defines the standard deviation (σ) of the Gaussian prior over lag:
%
%       lag ~ N(0, delayPriorSD^2)
%
%   This prior controls how much the model prefers lags near zero.
%
% Guidelines:
%
%   Scenario                      | Units     | Recommended delayPriorSD
%  ----------------------------------------------------------------------
%   Strong stimulus (e.g. CO₂)    | seconds   | 5 to 10
%   Resting-state with GM probe   | seconds   | 10 to 20
%   WM or noisy regions           | seconds   | 20 to 30
%   Unknown timing or free lag    | seconds   | 30+ (flat/uninformative)
%
% Notes:
%   • The prior is evaluated over `opts.lags` (which is also in seconds)
%   • A narrow prior biases the model toward near-zero lag
%   • A broad prior allows discovery of longer vascular delays
% ----------------------------------------------------------

function maps = LagVB(mask, gmMask, wmMask, data, probe, nuis, opts)
global opts;
opts.VBdir = fullfile(opts.resultsdir,'vbLAG'); mkdir(opts.VBdir);

if ~isfield(opts,'plotPrior'), opts.plotPrior = true; end
if ~isfield(opts,'delayPriorSD_GM'), opts.delayPriorSD_GM = 5; end %units of seconds
if ~isfield(opts,'delayPriorSD_WM'), opts.delayPriorSD_WM = 15; end
if ~isfield(opts,'lags'),     opts.lags    = -5:0.25:60; end
if ~isfield(opts,'polyOrd'),  opts.polyOrd = 1;          end
if ~isfield(opts,'priorVar'), opts.priorVar = 0.5;       end
if ~isfield(opts,'priorVar_GM'), opts.priorVar_GM = 5;       end   %(expected β value)^2  = (%BOLD/mmHg)²
if ~isfield(opts,'priorVar_WM'), opts.priorVar_WM = 2.5;       end
if ~isfield(opts,'priorDelay_GM'), opts.priorDelay_GM = 5;       end   %(expected β value)^2  = (%BOLD/mmHg)²
if ~isfield(opts,'priorDelay_WM'), opts.priorDelay_WM = 15;       end

if ~isfield(opts,'smoothmap'), opts.smoothmap = 0;       end
if ~isfield(opts,'scaleInput'), opts.scaleInput = 0;       end
if ~isfield(opts,'verbose'),  opts.verbose = 1;          end
if ~isfield(opts,'cvr'),  opts.cvr = 1;          end %compute CVR maps - makes sense if your probe is a CO2 trace

TR = opts.TR;
%% ---------- dimensions & sanity ---------------------------------------
% Tissue masks for adaptive priors
gmMask = logical(gmMask); wmMask = logical(wmMask);
tissueLabels = zeros(size(mask));
tissueLabels(gmMask) = 1;  % GM = 1
tissueLabels(wmMask) = 2;  % WM = 2
% Assign GM prior to any voxel in WB mask not covered by GM or WM
wbVoxels = find(mask);

for i = 1:length(wbVoxels)
    idx = wbVoxels(i);
    if tissueLabels(idx) == 0
        tissueLabels(idx) = 1;  % fallback to GM
    end
end

mask   = mask~=0;
[X,Y,Z,T] = size(data);
probe = double(probe(:)');
assert(T == numel(probe), 'Regressor length ≠ 4‑D time dimension');
if nargin < 5 || isempty(nuis), nuis = zeros(T,0); end

%% ---------- nuisance regressors ---------------------------------------
poly = []; for p = 0:opts.polyOrd, poly = [poly, ((1:T)' - T/2).^p]; end %#ok<AGROW>
poly = detrend(poly,'constant');
Xn   = [poly, nuis];
P    = size(Xn,2);

%% ---------- voxel list -------------------------------------------------
[voxTS, coords] = grabTimeseries(double(data), mask);
Nvox = size(voxTS,1);
if opts.verbose
    fprintf('lagVB | vox=%d  lag grid=%d  nuis=%d  TR=%.2f s\n', Nvox, numel(opts.lags), P, TR);
end

%% ---------- build shifted CO₂ matrix ----------------------------------
L = numel(opts.lags);
Xc = zeros(T,L);
ft = fft(probe).';
freq = (0:T-1)'/T;
for k = 1:L
    lag = opts.lags(k);
    phase = exp(-1i*2*pi*lag/TR .* freq);
    Xc(:,k) = real(ifft(ft .* phase));
end
% Check scaling: if CO2 is in mmHg, skip z-score to preserve units
if opts.scaleInput
    Xc = detrend(zscore(Xc), 'constant');  % default: z-score
else
    Xc = detrend(Xc, 'constant');  % preserve original units (e.g., mmHg)
end

%% ---------- constants --------------------------------------------------
% Prior mean for beta
mu0_GM = max(meanTimeseries(data,gmMask));
mu0_WM =  max(meanTimeseries(data,wmMask));  % or lower, depending on SNR and expected amplitude

mu_GM = opts.priorDelay_GM;       % GM lag prior mean in seconds
mu_WM = opts.priorDelay_WM;      % WM prior expects ~10s delay

% Prior mean for beta (CVR)
delayPrior_GM = normpdf(opts.lags, mu_GM, opts.delayPriorSD_GM);
delayPrior_WM = normpdf(opts.lags, mu_WM, opts.delayPriorSD_WM);

XtX_n = Xn' * Xn;
invV0 = 1/opts.priorVar;

if isfield(opts,'priorVar_GM'), invV0_GM = 1 / opts.priorVar_GM; end
if isfield(opts,'priorVar_WM'), invV0_WM = 1 / opts.priorVar_WM; else invV0_WM = invV0_GM; end

if opts.plotPrior
    figure('Name','Delay Prior','Color','w');
    plot(opts.lags, delayPrior_GM, 'm', 'LineWidth', 2); hold on;
    plot(opts.lags, delayPrior_WM, 'c', 'LineWidth', 2);
    xlabel('Lag (s)'); ylabel('Prior Probability'); grid on;
    legend('GM prior','WM prior');
    title('Delay Priors with Shifted Means');
    saveas(gcf,fullfile(opts.VBdir,'tissue_delay_priors.fig'));
end


%% ---------- allocate outputs ------------------------------------------
vec = nan(Nvox,1,'single');
[bestIdx, bestBeta, bestR2, bestSig2, postLag, postVar] = deal(vec);

%% ---------- PARFOR voxel loop -----------------------------------------
queue = createParallelProgressBar(numel(coords));  %  <---
parfor v = 1:Nvox
    label = tissueLabels(coords(v));
    if label == 1
        invV0 = invV0_GM; delayPrior = delayPrior_GM;
        mu0 = mu0_GM;
    elseif label == 2
        mu0 = mu0_WM;
        invV0 = invV0_WM; delayPrior = delayPrior_WM;
    else
        mu0 = mu0_GM;
        invV0 = invV0_GM; delayPrior = delayPrior_GM;
    end
    y = voxTS(v,:)';
    if all(y==0) || any(isnan(y)), continue; end

    Fvals = nan(L,1);
    betaVals = nan(L,1);

    bestF = -Inf; best_i = 1;
    best_beta = 0; best_sig2 = 0;
    best_m = zeros(P+1,1,'double'); best_xk = zeros(T,1,'double');

    for k = 1:L
        xk = Xc(:,k);
        X  = [xk, Xn];
        xkTXn = xk' * Xn;
        XtX = [xk'*xk, xkTXn ; xkTXn', XtX_n];
        Lambda0 = diag([invV0, zeros(1,P)]);
        S = XtX + Lambda0;
        ridge = 1e-6 * trace(S)/(P+1);
        S = S + eye(P+1)*ridge;
        [cholS,flag] = chol(S,'lower');
        if flag == 0
            logdet = 2*sum(log(diag(cholS)));
        else
            logdet = sum(log(max(eig(S), eps)));
        end
        b = X' * y + Lambda0 * [mu0; zeros(P,1)];
        m = S \ b;
        res = y - X*m;
        sig2 = (res' * res) / (T - (P+1));
        F = -0.5*(T*log(2*pi*sig2) + logdet + T) + log(delayPrior(k));

        Fvals(k) = F;
        betaVals(k) = m(1);

        if F > bestF
            bestF = F; best_i = k;
            best_beta = m(1);
            best_sig2 = sig2;
            best_m = m; best_xk = xk;
        end
    end

    bestIdx(v)  = best_i;
    bestBeta(v) = best_beta;
    bestSig2(v) = best_sig2;

    Xbest = [best_xk, Xn];
    yhat  = Xbest * best_m;
    SS_tot= sum((y-mean(y)).^2);
    SS_res= sum((y - yhat).^2);
    bestR2(v) = 1 - SS_res/SS_tot;

    % --- softmax posterior mean and variance of lag
    Fvals = Fvals - max(Fvals);
    w = exp(Fvals); w = w / sum(w);
    lagVals = opts.lags(:);
    postLag(v) = sum(w .* lagVals);
    postVar(v) = sum(w .* (lagVals - postLag(v)).^2);

    pause(0.00001);
    send(queue, v);
end

%% ---------- assemble maps ---------------------------------------------
mapSize = [X Y Z];
maps.tissueType = nan(mapSize,'single');
maps.tissueType(coords) = tissueLabels(coords);

lagVec = nan(Nvox,1,'single');
valid  = ~isnan(bestIdx);
lagVec(valid) = opts.lags(bestIdx(valid))';
values = {lagVec, bestBeta, bestR2, bestSig2, postLag, postVar};
fields = {'delay','beta','R2','sigma2','postLag','postVar'};
for i = 1:numel(fields)
    vol = nan(mapSize,'single');
    vol(coords) = values{i};
    maps.(fields{i}) = vol;
end

if opts.smoothmap
    temp = opts.FWHM;
    opts.FWHM = [5 5 5];
    maps.delay = filterData(single(maps.delay), mask, mask, opts);
    maps.postLag = filterData(single(maps.postLag), mask, mask, opts);
    opts.FWHM = temp; clear temp
end

tmp = maps.delay(:); tmp(tmp>0) = 0; offset = abs(min(tmp));
maps.delay = mask.*(maps.delay + offset);
tmp = maps.postLag(:); tmp(tmp>0) = 0; offset = abs(min(tmp));
maps.postLag = mask.*(maps.postLag + offset);

savedir = opts.VBdir;
saveMap(cast(mask.*maps.delay,opts.mapDatatype), savedir, 'delay', opts.info.map, opts);
saveMap(cast(mask.*maps.beta,opts.mapDatatype), savedir, 'beta', opts.info.map, opts);
saveMap(cast(mask.*maps.R2,opts.mapDatatype), savedir, 'R2', opts.info.map, opts);
saveMap(cast(mask.*maps.sigma2,opts.mapDatatype), savedir, 'sigma', opts.info.map, opts);
saveMap(cast(mask.*maps.postLag,opts.mapDatatype), savedir, 'postLag', opts.info.map, opts);
saveMap(cast(mask.*maps.postVar,opts.mapDatatype), savedir, 'postVar', opts.info.map, opts);
saveMap(cast(mask.*maps.tissueType,opts.mapDatatype), savedir, 'tissueType', opts.info.map, opts);

%% ---------- cvr mapping ------------------------------------------

if opts.cvr
    fprintf('Computing standard CVR and associated statistics...\n');

    % Allocate outputs
    cvr     = nan(Nvox,1);
    R2_CVR  = nan(Nvox,1);
    T_CVR   = nan(Nvox,1);

    probe = double(probe(:));
    T = length(probe);

    parfor v = 1:Nvox
        y = voxTS(v,:)';
        if all(y == 0) || any(isnan(y)), continue; end

        % Standard regression
        X = [probe, ones(T,1)];
        b = X \ y;
        yhat = X * b;

        cvr(v) = b(1);  % β = %BOLD/mmHg

        % R²
        SSE = norm(y - yhat)^2;
        SST = norm(y - mean(y))^2;
        R2_CVR(v) = 1 - (SSE / SST);

        % T-statistic
        resid_std = std(y - yhat);
        T_CVR(v) = b(1) / resid_std;
    end

    % Write to maps
    maps.cvr    = nan(mapSize, 'single');
    maps.R2CVR  = nan(mapSize, 'single');
    maps.TCVR   = nan(mapSize, 'single');

    maps.cvr(coords)   = single(cvr);
    maps.R2CVR(coords) = single(R2_CVR);
    maps.TCVR(coords)  = single(T_CVR);

    % Save
    saveMap(cast(mask .* maps.cvr,    opts.mapDatatype), opts.VBdir, 'cvr',    opts.info.map, opts);
    saveMap(cast(mask .* maps.R2CVR,  opts.mapDatatype), opts.VBdir, 'R2CVR',  opts.info.map, opts);
    saveMap(cast(mask .* maps.TCVR,   opts.mapDatatype), opts.VBdir, 'Tstat_CVR',   opts.info.map, opts);

    %%
    fprintf('Generating lag-adjusted CVR maps...\n');

    dyn = T;

    % Preallocate
    cvr_lag     = nan(Nvox,1);
    cR2     = nan(Nvox,1);
    cTstat  = nan(Nvox,1);

    shifted_regr = NaN(Nvox, dyn );

    % Compute shift indices (in TR units)
    shift_idx = round(maps.delay(coords) / TR);

    % Build shifted regressors per voxel
    for ii = 1:Nvox
        s = shift_idx(ii);                    % integer shift in TRs
        corr_regr = NaN(1, T);                % initialize shifted regressor

        if s == 0
            corr_regr = probe;
        elseif s > 0
            % Shift right (positive lag)
            corr_regr((s+1):end) = probe(1:end-s);
            corr_regr(1:s) = probe(1);        % pad beginning with first value
        else
            % Shift left (negative lag)
            s_abs = abs(s);
            corr_regr(1:end-s_abs) = probe((s_abs+1):end);
            corr_regr((end-s_abs+1):end) = probe(end);  % pad end with last value
        end

        shifted_regr(ii,:) = corr_regr;
    end

    % Fit voxelwise regressions
    regr_coef = zeros(2, Nvox);  % intercept + slope
    parfor ii = 1:Nvox
        A = shifted_regr(ii,:);
        B = voxTS(ii,:);
        B(isnan(A)) = []; A(isnan(A)) = [];
        if numel(A) < T/2, continue; end  % skip if mostly NaN

        X = [ones(length(A),1), A'];
        regr_coef(:,ii) = X \ B';
    end

    % Extract beta
    cvr_lag = regr_coef(2,:);
    cvr_lag(cvr_lag > 10 | cvr_lag < -10) = NaN;

    % Calculate R² and T-stats
    parfor ii = 1:Nvox
        A = shifted_regr(ii,:);
        B = voxTS(ii,:);
        B(isnan(A)) = []; A(isnan(A)) = [];
        if numel(A) < T/2, continue; end

        beta  = regr_coef(2,ii);
        alpha = regr_coef(1,ii);
        Y = alpha + beta * A';

        SSE = norm(B - Y')^2;
        SST = norm(B - mean(B))^2;

        cR2(ii)    = 1 - (SSE / SST);
        resid_std = std(B - Y');
        cTstat(ii)= beta / resid_std;
    end

    maps.cvrLagCorr = zeros([1, numel(mask)]);
    maps.R2LagCorr  = zeros([1, numel(mask)]);
    maps.TLagCorr   = zeros([1, numel(mask)]);

    maps.cvrLagCorr(coords) = single(cvr_lag);
    maps.R2LagCorr(coords)  = single(cR2);
    maps.TLagCorr(coords)   = single(cTstat);

    maps.cvrLagCorr = reshape(maps.cvrLagCorr, mapSize);
    maps.R2LagCorr = reshape(maps.R2LagCorr, mapSize);
    maps.TLagCorr = reshape(maps.TLagCorr, mapSize);

    % Save
    saveMap(mask .* maps.cvrLagCorr, opts.VBdir, 'cvrLagCorr', opts.info.map, opts);
    saveMap(mask .* maps.R2LagCorr, opts.VBdir, 'R2LagCorr', opts.info.map, opts);
    saveMap(mask .* maps.TLagCorr, opts.VBdir, 'Tstat_LagCorr',  opts.info.map, opts);
end
disp('finished running bayesian analysis')
disp('...saving maps in .mat file' )
save(fullfile(opts.VBdir,'result_lagVB_maps.mat'), 'maps');
end
