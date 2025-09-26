function [HRFidx, HRF_probe, HRFmeta] = convHRF(probe, opts)
% convHRF5 — FFT-style HRF convolution (matches convHRF2 mechanics)
%   • model: 'gamma' (default) or 'exp'
%   • use_onset: true/false  →  4 combinations total
%
% OUTPUTS
%   HRF_probe : (N × T) bank of regressors; row 1 = normalized probe
%   HRFidx    : (N × 3) indices = [onset_idx, disp_idx, under_idx]
%               onset_idx=0 denotes the zero-onset branch (no onset)
%               for 'exp', under_idx is set to 1 (N/A)
%   HRFmeta   : struct with grids and index meanings
%
% Notes
%   • Convolution is done in the frequency domain exactly like your convHRF2:
%       freq = fftshift(fft(kernel / trapz(kernel)));
%       tmp  = real(ifft(ifftshift( freq(1,end:-1:1) .* spec )));
%       row  = (tmp / trapz(tmp)); row = row(1,end:-1:1);
%   • Row 1 is the (area-normalized) padded probe, like your code.
%   • Padding/unpadding follow your pad/padfront behavior.
%   • Optional LD removal via licols if opts.removeLD is true.

warning('off');
global opts;

% ---------- defaults ----------
if ~isfield(opts,'model')      || isempty(opts.model),      opts.model = 'gamma'; end  % 'gamma'|'exp'
if ~isfield(opts,'use_onset')  || isempty(opts.use_onset),  opts.use_onset = true; end
if ~isfield(opts,'onset')      || isempty(opts.onset),      opts.onset = 1; end        % gamma: alpha1; exp: Δ (samples)
if ~isfield(opts,'disp')       || isempty(opts.disp),       opts.disp  = 1:0.5:30; end % gamma: beta1;  exp: τ (samples)
if ~isfield(opts,'under')      || isempty(opts.under),      opts.under = 1; end        % gamma only (beta2)
if ~isfield(opts,'rratio')     || isempty(opts.rratio),     opts.rratio = 1000; end
if ~isfield(opts,'pad')        || isempty(opts.pad),        opts.pad = true; end
if ~isfield(opts,'padfront')   || isempty(opts.padfront),   opts.padfront = true; end
if ~isfield(opts,'removeLD')   || isempty(opts.removeLD),   opts.removeLD = false; end
if ~isfield(opts,'tol')        || isempty(opts.tol),        opts.tol = 1e-10; end
if ~isfield(opts,'verbose')    || isempty(opts.verbose),    opts.verbose = false; end
if ~isfield(opts,'rescale01')  || isempty(opts.rescale01),  opts.rescale01 = true; end

% guard against empty grids
if isempty(opts.disp),   opts.disp   = 1; end
if isempty(opts.under),  opts.under  = 1; end
if isempty(opts.onset),  opts.onset  = 1; end

% ---------- helpers ----------
rowify    = @(x) x(:).';                         % force row
area_norm = @(x) x ./ max(trapz(x), eps);
gamma_pdf = @(tt, a, b) ((tt.^(a-1)) .* (b.^(-a)) .* exp(-tt./b)) ./ gamma(a);

% ---------- pad probe exactly like convHRF2 ----------
probe = probe(:);
input_probe = rescale(probe);                     % your initial rescale
T0 = numel(input_probe);
if opts.pad, pad = 2^nextpow2(T0); else, pad = 0; end

if opts.padfront && opts.pad
    probeP = zeros(T0+pad,1);
    probeP(1:pad)        = input_probe(1);
    probeP(pad+1:end)    = input_probe;
else
    probeP = zeros(T0+2*pad,1);
    probeP(1:pad)              = input_probe(1);
    probeP(pad+1:end-pad)      = input_probe;
    probeP(end-pad+1:end)      = input_probe(end);
end

T      = numel(probeP);
t      = (0:T-1);
rr     = max(opts.rratio,1);
rowset = {};     % collect rows here (as row-vectors)
idxset = [];     % [onset_idx, disp_idx, under_idx]

% ---------- probe spectrum (as in your code) ----------
newprobe_freq = fftshift(fft(probeP) / max(trapz(probeP), eps));

% ---------- Row 1: normalized probe (area), as in convHRF2 ----------
rowset{end+1} = rowify(probeP / max(trapz(probeP), eps));
idxset(end+1,:) = [1,1,1];   % your first-row tag

model = lower(string(opts.model));
useOn = logical(opts.use_onset);

switch model
% ============================== GAMMA =====================================
case "gamma"
    % Onset kernels: alpha1=onset, beta1=1; undershoot alpha2=1, beta2=0.5 (scaled 1/rratio)
    onset_list = opts.onset(:).';
    if ~useOn
        onset_list = onset_list(1);   % fixed onset; zero-onset handled below
    end

    % Dispersion family: alpha1=1, beta1=disp; undershoot alpha2=1, beta2=under (scaled 1/rratio)
    disp_list  = opts.disp(:).';
    under_list = opts.under(:).';

    % ----- build onset bank (time-domain, then normalize) -----
    HRFon = zeros(numel(onset_list), T);
    for ii = 1:numel(onset_list)
        a1 = onset_list(ii); b1 = 1;
        a2 = 1;              b2 = 0.5;
        h1 = gamma_pdf(t, a1, b1);
        h2 = gamma_pdf(t, a2, b2) / rr;
        hon = h1 - h2;
        HRFon(ii,:) = hon / max(trapz(hon), eps);
    end

    % ----- build dispersion bank -----
    HRFdsp = zeros(numel(disp_list)*numel(under_list), T);
    dsp_map = zeros(size(HRFdsp,1),2); % [disp_idx, under_idx]
    mm = 0;
    for jj = 1:numel(disp_list)
        for kk = 1:numel(under_list)
            mm = mm + 1;
            a1 = 1; b1 = disp_list(jj);
            a2 = 1; b2 = under_list(kk);
            h1 = gamma_pdf(t, a1, b1);
            h2 = gamma_pdf(t, a2, b2) / rr;
            h  = h1 - h2;
            HRFdsp(mm,:) = h / max(trapz(h), eps);
            dsp_map(mm,:) = [jj, kk];
        end
    end
    pp = size(HRFdsp,1);

    % ----- onset branch: probe ⊗ onset (FFT style), then ⊗ dispersion -----
    if useOn
        for ii = 1:size(HRFon,1)
            HRFon_freq = fftshift(fft(HRFon(ii,:) / max(trapz(HRFon(ii,:)),eps)));
            tmp_on     = real(ifft(ifftshift( HRFon_freq(1,end:-1:1) .* newprobe_freq' )));
            onset_sig  = tmp_on / max(trapz(tmp_on), eps);
            onfreq     = fftshift(fft(onset_sig) / max(trapz(probeP), eps));

            for jj = 1:pp
                HRF_freq  = fftshift(fft(HRFdsp(jj,:) / max(trapz(HRFdsp(jj,:)),eps)));
                tmp       = real(ifft(ifftshift( HRF_freq(1,end:-1:1) .* onfreq )));
                tmp2      = tmp / max(trapz(tmp), eps);
                row       = tmp2(1,end:-1:1);            % your time flip
                rowset{end+1} = rowify(row);
                idxset(end+1,:) = [ii, dsp_map(jj,1), dsp_map(jj,2)];
            end
        end
    end

    % ----- zero-onset branch: probe ⊗ dispersion only (exactly like your code) -----
    for jj = 1:pp
        HRF_freq = fftshift(fft(HRFdsp(jj,:) / max(trapz(HRFdsp(jj,:)),eps)));
        tmp      = real(ifft(ifftshift( HRF_freq(1,end:-1:1) .* newprobe_freq' )));
        tmp2     = tmp / max(trapz(tmp), eps);
        row      = tmp2(1,end:-1:1);
        rowset{end+1} = rowify(row);
        idxset(end+1,:) = [0, dsp_map(jj,1), dsp_map(jj,2)];
    end

% =============================== EXP ======================================
case "exp"
    % τ grid from opts.disp; optional onset delays Δ from opts.onset
    tau_list = opts.disp(:).';
    if useOn, onsetD = opts.onset(:).'; else, onsetD = 0; end

    % causal exponential kernels (in time), then FFT-style conv like above
    for io = 1:numel(onsetD)
        dlt = max(round(onsetD(io)), 0);
        for jj = 1:numel(tau_list)
            tau = max(tau_list(jj), eps);

            % make kernel in time, normalize by area
            k = zeros(1,T);
            if dlt < T
                k(dlt+1:end) = exp(-(0:T-dlt-1)/tau);
            end
            k = k / max(trapz(k), eps);

            kfreq = fftshift(fft(k));
            tmp   = real(ifft(ifftshift( kfreq(1,end:-1:1) .* newprobe_freq' )));
            tmp2  = tmp / max(trapz(tmp), eps);
            row   = tmp2(1,end:-1:1);

            rowset{end+1} = rowify(row);
            onset_idx = io; if ~useOn, onset_idx = 0; end
            idxset(end+1,:) = [onset_idx, jj, 1];   % under not used for EXP
        end
    end

otherwise
    error('opts.model must be ''gamma'' or ''exp''.');
end

% ---------- stack & unpad exactly like convHRF2 ----------
HRF_full = cell2mat(rowset.');   % all rows have identical length & shape
if opts.padfront && opts.pad
    HRF_probe = HRF_full(:, pad+1:end);
else
    HRF_probe = HRF_full(:, pad+1:end - double(opts.pad)*pad);
end
HRFidx = idxset;

% ---------- optional LD removal ----------
if opts.removeLD && size(HRF_probe,1) > 1
    [Xsub,keep] = licols(HRF_probe', opts.tol);  % your helper
    HRF_probe   = Xsub';
    HRFidx      = HRFidx(keep,:);
end

% ---------- final per-row rescale (your style) ----------
if opts.rescale01
    for ii=1:size(HRF_probe,1)
        HRF_probe(ii,:) = rescale(HRF_probe(ii,:), -1, 1);
    end
end

% ---------- meta ----------
HRFmeta = struct();
HRFmeta.model         = char(model);
HRFmeta.use_onset     = logical(useOn);
HRFmeta.onset_values  = opts.onset(:).';
HRFmeta.disp_values   = opts.disp(:).';
HRFmeta.under_values  = strcmpi(HRFmeta.model,'exp') * 0 + ~strcmpi(HRFmeta.model,'exp') * opts.under(:).';
HRFmeta.index_columns = struct('onset_idx',1,'disp_idx',2,'under_idx',3);
HRFmeta.notes         = 'Row 1 = normalized padded probe. onset_idx=0 denotes zero-onset branch. EXP ignores undershoot.';

% ---------- quick preview ----------
if opts.verbose
    figure('Units','pixels','Position',[220,520,720,190]);
    hold on; cm = parula(size(HRF_probe,1));
    for ii=1:size(HRF_probe,1), plot(HRF_probe(ii,:), 'Color', cm(ii,:)); end
    plot(HRF_probe(1,:), 'k', 'LineWidth', 2);
    title(sprintf('convHRF5: %s, onset=%d, N=%d', HRFmeta.model, HRFmeta.use_onset, size(HRF_probe,1)));
    xlim([1 size(HRF_probe,2)]); box off; hold off
end
end
