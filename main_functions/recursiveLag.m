% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% This implementation was developed by Stefan Rademakers,
% stefan-rademakers@outlook.com
% Sections of this code were contributed by Allen A. Champagne, a.champagne@queensu.ca
% The recursive lag approach is based on the original code shared by Dr.
% Toshiko Aso:
% https://github.com/aso-toshihiko/BOLDLagMapping_Deperfusioning
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>

% --- Documentation ---

% recursiveLag.m
%
% PURPOSE
%   Recursive estimation of hemodynamic/BOLD lag structure across the brain
%   (adapted from Aso et al., 2017). The algorithm starts from a seed
%   regressor and iteratively propagates lag information, updating the
%   regressor at each step with voxels that show strong local
%   cross‑correlation. Yseful for resting state lag mapping.
%
% USAGE
%   maps = recursiveLag(mask, data, probe, opts)
%
% INPUTS
%   mask      – 3‑D logical (X×Y×Z) brain mask.
%   data    – 4‑D numeric (X×Y×Z×T) fMRI/BOLD data.
%   probe – (T×1) reference timeseries.
%   opts      – Struct of optional parameters (see below).
%
% OUTPUT
%   maps.recursiveLag_map : Lag (seconds) for each voxel
%   maps.recursiveR_map   : Peak cross‑correlation value per voxel
%
% -------------------------------------------------------------------------
% opts FIELDS (defaults)
%   TR             (req)  – Repetition time (s)
%   interp_factor  (2)    – Temporal up‑sampling factor
%   lim            (3)    – Local ±lag search (samples after interp)
%   lim_s          (1)    – Stored subset of local lags
%   lowlag         (‑5)   – Minimum global lag (samples @ native TR)
%   highlag        (30)   – Maximum global lag
%   THR            (0.3)  – |r| threshold for accepting voxels
%   uni            (0)    – If 1, ignore negative correlations
%   prewhite       (0)    – Pre‑whiten time‑series before xcorr
%   rescale_probe  (1)    – Rescale evolving regressor to [0 1]
%   showwaitbar    (1)    – Display GUI progress bar
%   resultsdir     (pwd)  – Output folder for NIfTI maps
%

% ---------------------------------------------------------------
function maps = recursiveLag(mask, data, probe, opts)
global opts
%% Logical mask
mask = logical(mask);

%% ---------------- Parameter defaults -------------------------
if ~isfield(opts,'plot'),             opts.plot           = 1;  end
if ~isfield(opts,'prewhite'),         opts.prewhite       = 1;  end
if ~isfield(opts,'interp_factor'),    opts.interp_factor  = 1;  end
if ~isfield(opts,'uni'),              opts.uni            = 0;  end
if ~isfield(opts,'rescale_probe'),    opts.rescale_probe  = 0;  end
if ~isfield(opts,'norm_regr'),        opts.norm_regr      = 0;  end
if ~isfield(opts,'lim'),              opts.lim            = 3;  end
if ~isfield(opts,'lim_s'),            opts.lim_s          = 3;  end
if ~isfield(opts,'THR'),              opts.THR            = 0.3; end
if ~isfield(opts,'resultsdir'),       opts.resultsdir     = pwd;end
if ~isfield(opts,'lowlag'),           opts.lowlag         = -8; end
if ~isfield(opts,'highlag'),          opts.highlag        = 8; end
if ~isfield(opts,'TR'),              error('opts.TR (seconds) must be supplied');end
if ~isfield(opts,'showwaitbar'),     opts.showwaitbar    = 1;  end
if ~isfield(opts,'smoothmap'),       opts.smoothmap    = 1;  end
if ~isfield(opts,'circular'), opts.circular = 1; end   % true = use circshift

%% Derived parameters
opts.adjlowlag  = opts.lowlag  * opts.interp_factor;
opts.adjhighlag = opts.highlag * opts.interp_factor;
opts.recursivedir = fullfile(opts.resultsdir,'recursiveLag');
if ~exist(opts.recursivedir,'dir'), mkdir(opts.recursivedir); end

% ---------------- Data prep ------------------------
[xx,yy,zz,dyn] = size(data);
probe = probe(:);

[orig_ts, coordinates] = grabTimeseries(data, mask);

% Pre‑whiten if requested
if opts.prewhite
    pw_ts = orig_ts';
    parfor v = 1:numel(coordinates)
        pw_ts(:,v) = prewhiten(pw_ts(:,v));
    end
    pw_ts = pw_ts';
    probe = prewhiten(probe);
end

% Temporal interpolation
probe = interp(probe, opts.interp_factor);
wb_ts = zeros(numel(coordinates), opts.interp_factor*dyn, 'single');
if opts.prewhite
    parfor v = 1:numel(coordinates)
        wb_ts(v,:) = interp(pw_ts(v,:), opts.interp_factor);
    end
else
    parfor v = 1:numel(coordinates)
        wb_ts(v,:) = interp(orig_ts(v,:), opts.interp_factor);
    end
end

% Store evolving regressors for plotting
if opts.plot
    maxRegs     = 1 + 2*max(opts.adjhighlag-1, abs(opts.adjlowlag)-1);
    evolvingReg = cell(maxRegs,1); evolvingReg{1} = probe; regCounter = 1;
end

% Lag window setup
lag_window   = -opts.lim:opts.lim;    lag0LW = opts.lim+1;
keep_idx     = lag0LW + (-opts.lim_s:opts.lim_s);
Nlags        = opts.adjhighlag - opts.adjlowlag + 1;
col0         = 1 - opts.adjlowlag;   % global column for lag 0
r_map   = zeros(xx*yy*zz, Nlags,'single');
lag_map =  NaN (xx*yy*zz, Nlags,'single');

% Initial X‑corr
C0 = cell2mat(cellfun(@(a) xcorr(a,probe,opts.lim,'coeff'), mat2cell(wb_ts,ones(size(wb_ts,1),1),size(wb_ts,2)), 'uni',0));
if opts.uni, [R0,peak] = max(C0,[],2); else, [R0,peak] = max(abs(C0),[],2); end
R0(R0<opts.THR)=0;

for k = 1:numel(coordinates)
    if ismember(peak(k),keep_idx)
        r_map(coordinates(k),col0)   = R0(k);
        lag_map(coordinates(k),col0) = lag_window(peak(k));
    end
end

Uidx = peak; Lidx = peak; Rpos=R0; Rneg=R0;

% -------------------------------------------------------------------------
% Recursive propagation ----------------------------------------------------
% -------------------------------------------------------------------------
steps = max(opts.adjhighlag-1, abs(opts.adjlowlag)-1);
if opts.showwaitbar, hWB = waitbar(0, 'Recursive lag mapping...'); end

for p = 1:steps
    %% ---------------- Positive branch -----------------------------------
    if p <= opts.adjhighlag-1
        vox = (Uidx == lag0LW + 1) & (Rpos > opts.THR);
        if any(vox)
            % ─── POSITIVE branch regressor ─────────────────────────────────────────
            reg = mean(wb_ts(vox,:), 1);
            if opts.rescale_probe, reg = rescale(reg); end

            if opts.circular
                % Classic circular shift (no NaNs)
                Ureg = circshift(reg,  [0  p]);           % delay by +p samples
            else
                % NaN-pad then inpaint gaps
                Ureg          = nan(size(probe), 'like', probe);
                Ureg(1+p:end) = reg(1:end-p);
                Ureg          = inpaint_nans(Ureg);       % or fillgaps / fillmissing
            end

            % Diagnostics
            if opts.plot
                regCounter              = regCounter + 1;
                evolvingReg{regCounter} = Ureg;   % ← store the filled regressor
            end

            % Cross-correlation (unchanged)
            C = cell2mat(cellfun(@(a) xcorr(a, Ureg, opts.lim, 'coeff'), ...
                mat2cell(wb_ts, ones(size(wb_ts,1),1), size(wb_ts,2)), ...
                'UniformOutput', false));
            if opts.uni,  [Rpos, Uidx] = max(C, [], 2);
            else,         [Rpos, Uidx] = max(abs(C), [], 2); end
            Rpos(Rpos < opts.THR) = 0;

            col = col0 + p;   % global column index
            for k = 1:numel(coordinates)
                if ismember(Uidx(k), keep_idx)
                    r_map(coordinates(k), col)   = Rpos(k);
                    lag_map(coordinates(k), col) = lag_window(Uidx(k)) + p;
                end
            end
        end
    end

    %% ---------------- Negative branch -----------------------------------
    if p <= abs(opts.adjlowlag) - 1
        vox = (Lidx == lag0LW - 1) & (Rneg > opts.THR);
        if any(vox)
            % ─── NEGATIVE branch regressor ────────────────────────────────────────
            reg = mean(wb_ts(vox,:), 1);

            if opts.circular
                Lreg = circshift(reg,  [0 -p]);           % advance by –p samples
            else
                Lreg          = nan(size(probe), 'like', probe);
                Lreg(1:end-p) = reg(1+p:end);
                Lreg          = inpaint_nans(Lreg);
            end

            C = cell2mat(cellfun(@(a) xcorr(a, Lreg, opts.lim, 'coeff'), ...
                mat2cell(wb_ts, ones(size(wb_ts,1),1), size(wb_ts,2)), ...
                'UniformOutput', false));
            if opts.uni,  [Rneg, Lidx] = max(C, [], 2);
            else,         [Rneg, Lidx] = max(abs(C), [], 2); end
            Rneg(Rneg < opts.THR) = 0;

            col = col0 - p;
            for k = 1:numel(coordinates)
                if ismember(Lidx(k), keep_idx)
                    r_map(coordinates(k), col)   = Rneg(k);
                    lag_map(coordinates(k), col) = lag_window(Lidx(k)) - p;
                end
            end
        end
    end

    % Update waitbar
    if opts.showwaitbar && ishandle(hWB)
        waitbar(p/steps, hWB);
    end
end

if opts.showwaitbar && exist('hWB', 'var') && ishandle(hWB)
    close(hWB);
end

% -------------------------------------------------------------------------
% Extract max‑r lag per voxel ---------------------------------------------
% -------------------------------------------------------------------------
[maxR, bestCol] = max(r_map, [], 2);
lag_best        = NaN(xx*yy*zz, 1, 'single');
for k = 1:numel(coordinates)
    lag_best(coordinates(k)) = lag_map(coordinates(k), bestCol(coordinates(k)));
end

lag_sec = opts.TR * lag_best / opts.interp_factor;

% -------------------------------------------------------------------------
% Diagnostic plot of evolving regressors ----------------------------------
% -------------------------------------------------------------------------
if opts.plot && regCounter >= 2
    evolvingReg = evolvingReg(1:regCounter);
    taxis       = (0:numel(evolvingReg{1})-1) * opts.TR / opts.interp_factor;

    figure('Name', 'Evolving regressors', 'Color', 'w'); hold on
    cmap = jet(numel(evolvingReg));
    for ii = 1:numel(evolvingReg)
        plot(taxis, zscore(double(evolvingReg{ii})), 'Color', cmap(ii,:));
    end
    xlabel('Time (s)'); ylabel('z‑scored amplitude'); grid on
    title('Average time‑course of each updated regressor');
    colormap(cmap); cb = colorbar; cb.Label.String = 'Iteration (early → late)';
end

% -------------------------------------------------------------------------
% Prepare output maps ------------------------------------------------------
% -------------------------------------------------------------------------
lag_img = mask .* reshape(lag_sec, [xx yy zz]);
R_img   = mask .* reshape(maxR  , [xx yy zz]);

maps.recursiveLag_map = lag_img;
maps.recursiveR_map   = R_img;
savedir = opts.recursivedir;

% Optional spatial smoothing of lag map
if opts.smoothmap
    fwhm_tmp  = opts.FWHM;        % backup
    opts.FWHM = [3 3 3];
    maps.recursiveLag_map_smoothed = filterData(single(maps.recursiveLag_map), mask, mask, opts);
    opts.FWHM = fwhm_tmp; clear fwhm_tmp
    saveMap(cast(maps.recursiveLag_map_smoothed, opts.mapDatatype), savedir, 'recursive_lag_map_smoothed', opts.info.map, opts);
end

% -------------------------------------------------------------------------
% Save maps ---------------------------------------------------------------
% -------------------------------------------------------------------------
if ~exist(savedir, 'dir'), mkdir(savedir); end
saveMap(cast(lag_img, opts.mapDatatype), savedir, 'recursive_lag_map', opts.info.map, opts);
saveMap(cast(R_img,   opts.mapDatatype), savedir, 'recursive_r_map',   opts.info.map, opts);

end