% =====================================================================
% lagVB.m  –  Bayesian (grid) lag mapping for CVR
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
%   maps = lagVB(mask, data, probe, TR, nuis, opts)
% ---------------------------------------------------------------------
% Outputs (single‑precision maps)
%   maps.delay    – best lag [s]
%   maps.beta     – β at best lag
%   maps.R2       – R² of full model at best lag
%   maps.sigma2   – residual variance
%   maps.postLag  – posterior mean lag [s]
%   maps.postVar  – posterior variance of lag [s²]
% ---------------------------------------------------------------------
function maps = LagVB(mask, gmMask, wmMask, data, probe, nuis, opts)
global opts;
opts.VBdir = fullfile(opts.resultsdir,'vbLAG'); mkdir(opts.VBdir);

%% ---------- defaults ---------------------------------------------------
% Prior configuration
% - opts.priorVar controls the prior variance on β (CVR amplitude). A large value (e.g. 2000) weakly constrains β, useful in GM.
% - opts.delayPriorSD sets the standard deviation (in seconds) of the Gaussian delay prior. Smaller values bias toward near-zero lags.
% - For tissue-adaptive priors, you can provide:
%     opts.priorVar_GM: variance for gray matter
%     opts.priorVar_WM: variance for white matter
%     opts.delayPriorSD_GM: delay SD for gray matter
%     opts.delayPriorSD_WM: delay SD for white matter

%priorVar = 1000 → very weak constraint, wide prior, good for flexible GM fits.
%priorVar = 100 → moderate constraint, good balance between prior and data.
%priorVar = 10 → tight constraint, useful in low SNR voxels or WM to dampen β inflation.

% The function will fall back to GM priors for voxels outside GM or WM masks.

if ~isfield(opts,'plotPrior'), opts.plotPrior = false; end
if ~isfield(opts,'delayPriorSD'), opts.delayPriorSD = 10; end
if ~isfield(opts,'delayPriorSD_GM'), opts.delayPriorSD_GM = 5; end
if ~isfield(opts,'delayPriorSD_WM'), opts.delayPriorSD_WM = 30; end
if ~isfield(opts,'lags'),     opts.lags    = -5:0.25:60; end
if ~isfield(opts,'polyOrd'),  opts.polyOrd = 1;          end
if ~isfield(opts,'priorVar'), opts.priorVar = 800;       end
if ~isfield(opts,'priorVar_GM'), opts.priorVar_GM = 500;       end
if ~isfield(opts,'priorVar_WM'), opts.priorVar_WM = 50;       end
if ~isfield(opts,'smoothmap'), opts.smoothmap = 1;       end

if ~isfield(opts,'verbose'),  opts.verbose = 1;          end

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
if isfield(opts,'scaleInput') && opts.scaleInput == false
    Xc = detrend(Xc, 'constant');  % preserve original units (e.g., mmHg)
else
    Xc = detrend(zscore(Xc), 'constant');  % default: z-score
end

%% ---------- constants --------------------------------------------------
mu0 = 0.5;  % Prior mean for beta
mu0_GM = 0.5;
mu0_WM = 0.2;  % or lower, depending on SNR and expected amplitude

% Prior mean for beta (CVR)
delayPrior_GM = normpdf(opts.lags, 0, opts.delayPriorSD);  % fallback default
if isfield(opts,'delayPriorSD_GM'), delayPrior_GM = normpdf(opts.lags, 0, opts.delayPriorSD_GM); end
if isfield(opts,'delayPriorSD_WM'), delayPrior_WM = normpdf(opts.lags, 0, opts.delayPriorSD_WM); else delayPrior_WM = delayPrior_GM; end
invV0_GM = 1 / opts.priorVar;
if isfield(opts,'priorVar_GM'), invV0_GM = 1 / opts.priorVar_GM; end
if isfield(opts,'priorVar_WM'), invV0_WM = 1 / opts.priorVar_WM; else invV0_WM = invV0_GM; end

delayPrior = normpdf(opts.lags, 0, opts.delayPriorSD);  % Prior over delays

if opts.plotPrior
    figure('Name','Delay Prior','Color','w');
    plot(opts.lags, delayPrior, 'k', 'LineWidth', 2);
    xlabel('Lag (s)'); ylabel('Prior Probability'); grid on;
    title(sprintf('Delay Prior: N(0, %g^2)', opts.delayPriorSD));
end
XtX_n = Xn' * Xn;
invV0 = 1/opts.priorVar;

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
opts.FWHM = [3 3 3];
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

end
