% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <convHRF: convolved an input signal with a exponential response function >
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
%
% ************************************************************************
function [HRFidx, HRF_probe] = convHRF2(probe, opts)
% convHRF2
% Convolve input probe with a bank of double-gamma HRFs varying onset & dispersion.
% OUTPUTS:
%   HRFidx   : [onset_idx, disp_idx, 1] per row of HRF_probe. onset_idx=0 => zero-onset branch.
%   HRF_probe: rows are the normalized probe convolved with each HRF variant (row 1 = normalized probe)
global opts

% ----------------------- defaults (non-destructive) -----------------------
opts = setDefault(opts, 'verbose',          false);      % plot or not
opts = setDefault(opts, 'rratio',           1000);       % large => small undershoot (1/rratio)
opts = setDefault(opts, 'onset',            1);          % gamma shape(s) for main & undershoot 
opts = setDefault(opts, 'disp',             1:0.5:30);   % dispersion (samples), kept compatible
opts = setDefault(opts, 'pad',              true);
opts = setDefault(opts, 'padfront',         true);
opts = setDefault(opts, 'detrendHRF',       false);      % keep your flag
opts = setDefault(opts, 'TR',               1.0);        % seconds (used only for plotting/area calc if desired)
opts = setDefault(opts, 'norm_mode',        'area');     % 'area'|'l2'|'none' (final row-wise normalization)
opts = setDefault(opts, 'tol',              1e-10);      % rank reduction tolerance
% Optional: allow exponential dispersion to be in seconds, if user supplies this:
opts = setDefault(opts, 'disp_in_seconds',  false);      % if true, interprets opts.disp as seconds (uses opts.TR)
% -------------------------------------------------------------------------

probe = probe(:);
N     = numel(probe);

% Initial rescale (as in your code)
input_probe = rescale(probe);

% ------------------------------ padding -----------------------------------
if opts.pad
    pad = 2^nextpow2(N);
else
    pad = 0;
end

if opts.padfront && opts.pad
    probePadded = zeros(N + pad, 1);
    probePadded(1:pad)        = input_probe(1);
    probePadded(pad+1:end)    = input_probe;
else
    probePadded = zeros(N + 2*pad, 1);
    probePadded(1:pad)              = input_probe(1);
    probePadded(pad+1:end-pad)      = input_probe;
    probePadded(end-pad+1:end)      = input_probe(end);
end

Np = numel(probePadded);

% Time vectors (samples and seconds)
t_samples = (0:Np-1).';                 % used to keep compatibility with your original dispersion usage
t_seconds = t_samples * opts.TR;

% ---------------------- build onset HRFs (double-gamma) -------------------
% Your original used alpha_1 = onset(ii), beta_1 = 1, alpha_2 = same alpha, beta_2 = 1,
% and undershoot scaled by 1/rratio.
onset_vals = opts.onset(:).';
nOn        = numel(onset_vals);
HRFon      = zeros(nOn, Np);

for ii = 1:nOn
    alpha_1 = onset_vals(ii);
    beta_1  = 1;
    alpha_2 = alpha_1;
    beta_2  = 1;

    % gamma lobes in *samples* to preserve original behavior
    h1 = gamma_pdf(t_samples, alpha_1, beta_1);
    h2 = gamma_pdf(t_samples, alpha_2, beta_2) / max(opts.rratio, eps);

    h_on = h1 - h2;
    h_on = h_on ./ max(trapz(h_on), eps);   % area-normalize
    HRFon(ii, :) = h_on.';
end

% ---------------------- build dispersion kernels --------------------------
% Original uses exp(-x/disp) with x in *samples*. Keep that default.
% If disp_in_seconds=true, interpret opts.disp as seconds.
if opts.disp_in_seconds
    tau = max(opts.disp(:), eps);                 % seconds
    HRFdsp = zeros(numel(tau), Np);
    for jj = 1:numel(tau)
        d = exp(-t_seconds / tau(jj));
        d = d ./ max(trapz(d), eps);
        HRFdsp(jj, :) = d.';
    end
else
    disp_vals = max(opts.disp(:), eps);           % samples
    HRFdsp = zeros(numel(disp_vals), Np);
    for jj = 1:numel(disp_vals)
        d = exp(-t_samples / disp_vals(jj));
        d = d ./ max(trapz(d), eps);
        HRFdsp(jj, :) = d.';
    end
end

nDs = size(HRFdsp,1);

% -------------------- allocate result matrices ----------------------------
% Variants = (nOn*nDs) + nDs + 1  : (onset⊗disp) + (probe⊗disp, zero-onset branch) + (probe itself)
nVar = (nOn * nDs) + nDs + 1;

HRF_probe_full = zeros(nVar, Np);
HRFidx_full    = zeros(nVar, 3);

% -------------------- row 1: probe itself (area-normalized) ---------------
k = 1;
HRFidx_full(k, :)    = [1, 1, 1];  % keep consistent with your earlier first set marker
probe_norm_area      = probePadded ./ max(trapz(probePadded), eps);
HRF_probe_full(k, :) = probe_norm_area.';
k = k + 1;

% ------------------ onset ⊗ dispersion, convolved with probe --------------
% We construct h_full = conv(onset, dispersion) (same length as probePadded),
% then y = conv(probePadded, h_full) with "same" cropping (no fftshift).
for io = 1:nOn
    h_on = HRFon(io, :).';
    for id = 1:nDs
        d      = HRFdsp(id, :).';
        h_full = convSame(h_on, d);                 % “same” length = Np

        y = convSame(probePadded, h_full);          % convolve with probe
        y = y ./ max(trapz(y), eps);                % match your normalization style

        HRF_probe_full(k, :) = y.';
        HRFidx_full(k, :)    = [io, id, 1];
        k = k + 1;
    end
end

% --------------- zero-onset branch: probe ⊗ dispersion only ---------------
for id = 1:nDs
    d  = HRFdsp(id, :).';
    y  = convSame(probePadded, d);
    y  = y ./ max(trapz(y), eps);

    HRF_probe_full(k, :) = y.';
    HRFidx_full(k, :)    = [0, id, 1];
    k = k + 1;
end

% ------------------------------ unpad -------------------------------------
if opts.padfront && opts.pad
    HRF_probe = HRF_probe_full(:, pad+1:end);
else
    HRF_probe = HRF_probe_full(:, pad+1:end-pad);
end
HRFidx = HRFidx_full;

% -------------------- rank reduction (QR, no licols) ----------------------
% Keep only linearly independent rows up to tolerance, preserving order.
% (Apply on the *rows* of HRF_probe.)
[keep] = rankReduceKeepIndices(HRF_probe, opts.tol);
HRF_probe = HRF_probe(keep, :);
HRFidx    = HRFidx(keep,    :);

% ------------------ optional detrend (if function exists) -----------------
if isfield(opts, 'detrendHRF') && opts.detrendHRF
    if exist('detrendHRF', 'file') == 2
        % Keep your original windowing fallback; adjust as needed.
        try
            HRF_probe = detrendHRF(HRF_probe, [min(50, size(HRF_probe,2)), max(size(HRF_probe,2)-20,1)]);
        catch
            warning('detrendHRF failed; skipping detrend.');
        end
    else
        warning('opts.detrendHRF=true but detrendHRF.m not found on path; skipping detrend.');
    end
end

% -------------------- final per-row normalization -------------------------
switch lower(string(opts.norm_mode))
    case "area"
        for ii = 1:size(HRF_probe,1)
            a = trapz(HRF_probe(ii, :));
            HRF_probe(ii, :) = HRF_probe(ii, :) ./ max(abs(a), eps);
        end
    case "l2"
        for ii = 1:size(HRF_probe,1)
            nrm = norm(HRF_probe(ii, :));
            HRF_probe(ii, :) = HRF_probe(ii, :) ./ max(nrm, eps);
        end
    otherwise
        % 'none' -> do nothing
end

% Finally rescale to [0,1] like your original
for ii = 1:size(HRF_probe,1)
    HRF_probe(ii, :) = rescale(HRF_probe(ii, :));
end

% ------------------------------ plot (optional) ---------------------------
if opts.verbose
    figure('Units','pixels','Position',[200, 500, 600, 160]);
    cmap = parula(size(HRF_probe,1));
    hold on;
    for ii = 1:size(HRF_probe,1)
        plot(HRF_probe(ii, :), 'Color', cmap(ii, :));
    end
    plot(HRF_probe(1, :), 'k', 'LineWidth', 2);  % emphasize row 1
    title('HRF-convolved probe'); xlim([1, size(HRF_probe,2)]); box off;

    if isfield(opts,'figdir') && ~isempty(opts.figdir)
        saveas(gcf, fullfile(opts.figdir, 'HRF_fnc.fig'));
    else
        saveas(gcf, fullfile(pwd, 'HRF_fnc.fig'));
    end
end

end % function convHRF2

% ============================== helpers ===================================

function y = convSame(x, h)
% Linear convolution with “same” length as x (centered crop)
    x = x(:); h = h(:);
    L  = numel(x) + numel(h) - 1;
    nfft = 2^nextpow2(L);
    yfull = ifft( fft(x, nfft) .* fft(h, nfft) );
    yfull = real(yfull(1:L));

    s = floor((numel(h)-1)/2);
    e = s + numel(x) - 1;
    y = yfull(s+1:e+1);
end

function g = gamma_pdf(t, k, theta)
% Explicit gamma PDF (no Statistics Toolbox), shape k>0, scale theta>0
% g(t) = t^(k-1) * exp(-t/theta) / (gamma(k) * theta^k), for t>=0; else 0
    t = max(t, 0);
    k = max(k, eps);
    theta = max(theta, eps);
    g = (t.^(k-1)) .* exp(-t./theta) ./ (gamma(k) .* (theta.^k));
    g(~isfinite(g)) = 0;
end

function s = setDefault(s, field, val)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = val;
    end
end

function keep = rankReduceKeepIndices(A, tol)
% Keep indices of linearly independent rows of A (row space),
% using QR on A' so we test column independence of A'.
% Returns indices in original order.
    if nargin < 2, tol = 1e-10; end
    [Q, R, E] = qr(A', 0);  %#ok<ASGLU>
    if isempty(R)
        keep = 1:size(A,1);
        return;
    end
    diagR = abs(diag(R));
    if isempty(diagR)
        keep = 1:size(A,1);
        return;
    end
    r = find(diagR > tol * max(1, diagR(1)), 1, 'last');
    if isempty(r), r = 1; end
    keep_unsorted = E(1:r);
    keep = sort(keep_unsorted); % preserve original order
end
