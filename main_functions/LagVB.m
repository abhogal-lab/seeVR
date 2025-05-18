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
%   maps = lagVB(mask, data4d, petco2, TR, nuis, opts)
% ---------------------------------------------------------------------
% Outputs (single‑precision maps)
%   maps.delay    – best lag [s]
%   maps.beta     – β at best lag
%   maps.R2       – R² of full model at best lag
%   maps.sigma2   – residual variance
%   maps.postLag  – posterior mean lag [s]
%   maps.postVar  – posterior variance of lag [s²]
% ---------------------------------------------------------------------
function maps = lagVB(mask, data4d, petco2, TR, nuis, opts)
global opts;
opts.VBdir = fullfile(opts.resultsdir,'vbLAG'); mkdir(opts.VBdir);
%% ---------- defaults ---------------------------------------------------
if ~isfield(opts,'lags'),     opts.lags    = -5:0.25:50; end
if ~isfield(opts,'polyOrd'),  opts.polyOrd = 1;          end
if ~isfield(opts,'priorVar'), opts.priorVar = 1e1;       end
if ~isfield(opts,'verbose'),  opts.verbose = 1;          end

%% ---------- dimensions & sanity ---------------------------------------
mask   = mask~=0;
[X,Y,Z,T] = size(data4d);
petco2 = double(petco2(:)');
assert(T == numel(petco2), 'Regressor length ≠ 4‑D time dimension');
if nargin < 5 || isempty(nuis), nuis = zeros(T,0); end

%% ---------- nuisance regressors ---------------------------------------
poly = []; for p = 0:opts.polyOrd, poly = [poly, ((1:T)' - T/2).^p]; end %#ok<AGROW>
poly = detrend(poly,'constant');
Xn   = [poly, nuis];
P    = size(Xn,2);

%% ---------- voxel list -------------------------------------------------
[voxTS, coords] = grabTimeseries(double(data4d), mask);
Nvox = size(voxTS,1);
if opts.verbose
    fprintf('lagVB | vox=%d  lag grid=%d  nuis=%d  TR=%.2f s\n', Nvox, numel(opts.lags), P, TR);
end

%% ---------- build shifted CO₂ matrix ----------------------------------
L = numel(opts.lags);
Xc = zeros(T,L);
ft = fft(petco2).';
freq = (0:T-1)'/T;
for k = 1:L
    lag = opts.lags(k);
    phase = exp(-1i*2*pi*lag/TR .* freq);
    Xc(:,k) = real(ifft(ft .* phase));
end
Xc = detrend(zscore(Xc), 'constant');

%% ---------- constants --------------------------------------------------
XtX_n = Xn' * Xn;
invV0 = 1/opts.priorVar;

%% ---------- allocate outputs ------------------------------------------
vec = nan(Nvox,1,'single');
[bestIdx, bestBeta, bestR2, bestSig2, postLag, postVar] = deal(vec);

%% ---------- PARFOR voxel loop -----------------------------------------
parfor v = 1:Nvox
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
        m = S \ (X' * y);
        res = y - X*m;
        sig2 = (res' * res) / (T - (P+1));
        F = -0.5*(T*log(2*pi*sig2) + logdet + T);

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
end

%% ---------- assemble maps ---------------------------------------------
mapSize = [X Y Z];
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

%shift delay to start at zero lag
maps.delay = maps.delay + abs(min(maps.delay));

if opts.niiwrite
    cd(opts.VBdir)
    niftiwrite(cast(mask.*maps.delay,opts.mapDatatype),'delay',opts.info.map);
    niftiwrite(cast(mask.*maps.beta,opts.mapDatatype),'beta',opts.info.map);
    niftiwrite(cast(mask.*maps.R2,opts.mapDatatype),'R2',opts.info.map);
    niftiwrite(cast(mask.*maps.sigma2,opts.mapDatatype),'sigma',opts.info.map);
    niftiwrite(cast(mask.*maps.postLag,opts.mapDatatype),'postLag',opts.info.map);
    niftiwrite(cast(mask.*maps.postVar,opts.mapDatatype),'postVar',opts.info.map);
else
    saveImageData(mask.*maps.delay, opts.headers.map, opts.VBdir, 'delay.nii.gz', datatype);
    saveImageData(mask.*maps.beta, opts.headers.map, opts.VBdir, 'beta.nii.gz', datatype);
    saveImageData(mask.*maps.R2, opts.headers.map, opts.VBdir, 'R2.nii.gz', datatype);
    saveImageData(mask.*maps.sigma2, opts.headers.map, opts.VBdir, 'sigma.nii.gz', datatype);
    saveImageData(mask.*maps.postLag, opts.headers.map, opts.VBdir, 'postLag.nii.gz', datatype);
    saveImageData(mask.*maps.postVar, opts.headers.map, opts.VBdir, 'postVar.nii.gz', datatype);
end

end