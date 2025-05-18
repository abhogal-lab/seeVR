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
%   cross‑correlation.
%
% USAGE
%   maps = recursiveLag(mask, data4D, seedProbe, opts)
%
% INPUTS
%   mask      – 3‑D logical (X×Y×Z) brain mask.
%   data4D    – 4‑D numeric (X×Y×Z×T) fMRI/BOLD data.
%   seedProbe – (T×1) reference timeseries.
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
if ~isfield(opts,'prewhite'),         opts.prewhite       = 0;  end
if ~isfield(opts,'interp_factor'),    opts.interp_factor  = 1;  end
if ~isfield(opts,'uni'),              opts.uni            = 0;  end
if ~isfield(opts,'rescale_probe'),    opts.rescale_probe  = 1;  end
if ~isfield(opts,'norm_regr'),        opts.norm_regr      = 0;  end
if ~isfield(opts,'lim'),              opts.lim            = 3;  end
if ~isfield(opts,'lim_s'),            opts.lim_s          = 1;  end
if ~isfield(opts,'THR'),              opts.THR            = 0.3;end
if ~isfield(opts,'comp'),             opts.comp           = 1;  end
if ~isfield(opts,'resultsdir'),       opts.resultsdir     = pwd;end
if ~isfield(opts,'lowlag'),           opts.lowlag         = -5; end
if ~isfield(opts,'highlag'),          opts.highlag        = 30; end
if ~isfield(opts,'TR'),              error('opts.TR (seconds) must be supplied');end
if ~isfield(opts,'showwaitbar'),     opts.showwaitbar    = 1;  end
if ~isfield(opts,'smoothmap'),       opts.smoothmap    = 1;  end


%% Derived parameters
opts.adjlowlag  = opts.lowlag  * opts.interp_factor;
opts.adjhighlag = opts.highlag * opts.interp_factor;
opts.recursivedir = fullfile(opts.resultsdir,'recursiveLag');
if ~exist(opts.recursivedir,'dir'), mkdir(opts.recursivedir); end

%% ---------------- Data prep ----------------------------------
[xx,yy,zz,dyn] = size(data);
probe = probe(:);

maps = struct();

[orig_voxel_ts, coordinates] = grabTimeseries(data, mask);

%% Optional pre‑whitening
if opts.prewhite
    pw_voxel_ts = orig_voxel_ts';
    parfor v = 1:size(coordinates,1)
        pw_voxel_ts(:,v) = prewhiten(pw_voxel_ts(:,v));
    end
    pw_voxel_ts = pw_voxel_ts';
    probe = prewhiten(probe);
end

%% Temporal interpolation
probe = interp(probe, opts.interp_factor);
wb_voxel_ts = zeros(length(coordinates), opts.interp_factor*dyn, 'single');
if opts.prewhite
    parfor v = 1:length(coordinates)
        wb_voxel_ts(v,:) = interp(pw_voxel_ts(v,:), opts.interp_factor);
    end
else
    parfor v = 1:length(coordinates)
        wb_voxel_ts(v,:) = interp(orig_voxel_ts(v,:), opts.interp_factor);
    end
end

%% ───────── OPTIONAL: capture & plot evolving regressors ─────────
if opts.plot
    % Pre-allocate a cell array large enough to hold all regressors
    maxRegs        = 1 + 2*max(opts.adjhighlag-1 , abs(opts.adjlowlag)-1);
    evolvingReg    = cell(maxRegs,1);
    regCounter     = 1;
    evolvingReg{regCounter} = probe;           % initial seed
end

%% Lag‑window definitions
lag_window   = -opts.lim:opts.lim;          % local search lags
lag_centerLW = opts.lim + 1;                % zero‑lag index in local window
store_range  = -opts.lim_s:opts.lim_s;
store_idx    = lag_centerLW + store_range;

nTotalLags   = opts.adjhighlag - opts.adjlowlag + 1;
centerCol    = 1 - opts.adjlowlag;          % global table column of lag 0

r_map   = zeros(xx*yy*zz, nTotalLags, 'single');
lag_map =  NaN(xx*yy*zz, nTotalLags, 'single');

%% ---------------- Initial correlation ------------------------
a2 = mat2cell(wb_voxel_ts, ones(size(wb_voxel_ts,1),1), size(wb_voxel_ts,2));
xcfun = @(a) xcorr(a, probe, opts.lim, 'coeff');
C0 = cell2mat(cellfun(xcfun, a2,'UniformOutput',false));
if opts.uni, [R0, idx0] = max(C0,[],2); else, [R0, idx0] = max(abs(C0),[],2); end
R0(R0<opts.THR) = 0;

lag_vals = lag_window;
for k = 1:length(coordinates)
    row = coordinates(k);
    li  = idx0(k);
    if ismember(li, store_idx)
        r_map(row, centerCol)   = R0(k);
        lag_map(row, centerCol) = lag_vals(li);
    end
end

Uidx = idx0; Lidx = idx0; Rpos = R0; Rneg = R0;

%% ---------------- Recursive propagation ----------------------
maxSteps = max(opts.adjhighlag-1 , abs(opts.adjlowlag)-1);
if opts.showwaitbar, hWB = waitbar(0,'Recursive lag mapping...'); end

for p = 1:maxSteps
    % -------- Positive branch --------
    if p <= opts.adjhighlag-1
        voxSel = (Uidx == lag_centerLW+1) & (Rpos > opts.THR);
        if any(voxSel)
            reg = mean(wb_voxel_ts(voxSel,:),1);
            if opts.rescale_probe, reg = rescale(reg); end
            Uregr = zeros(size(probe),'like',probe);
            Uregr(1+p:end) = reg(1:end-p);
            if opts.plot
                regCounter = regCounter + 1;
                evolvingReg{regCounter} = Uregr;
            end
            C = cell2mat(cellfun(@(a) xcorr(a,Uregr,opts.lim,'coeff'), a2,'UniformOutput',false));
            if opts.uni, [Rpos,Uidx] = max(C,[],2); else, [Rpos,Uidx] = max(abs(C),[],2); end
            Rpos(Rpos<opts.THR) = 0;

            col = centerCol + p;
            for k = 1:length(coordinates)
                row = coordinates(k); li = Uidx(k);
                if ismember(li, store_idx)
                    r_map(row,col)   = Rpos(k);
                    lag_map(row,col) = lag_vals(li) + p;
                end
            end
        end
    end

    % -------- Negative branch --------
    if p <= abs(opts.adjlowlag)-1
        voxSel = (Lidx == lag_centerLW-1) & (Rneg > opts.THR);
        if any(voxSel)
            reg = mean(wb_voxel_ts(voxSel,:),1);
            Lregr = zeros(size(probe),'like',probe);
            Lregr(1:end-p) = reg(1+p:end);
            if opts.plot
                regCounter = regCounter + 1;
                evolvingReg{regCounter} = Lregr;
            end
            C = cell2mat(cellfun(@(a) xcorr(a,Lregr,opts.lim,'coeff'), a2,'UniformOutput',false));
            if opts.uni, [Rneg,Lidx] = max(C,[],2); else, [Rneg,Lidx] = max(abs(C),[],2); end
            Rneg(Rneg<opts.THR) = 0;

            col = centerCol - p;
            for k = 1:length(coordinates)
                row = coordinates(k); li = Lidx(k);
                if ismember(li, store_idx)
                    r_map(row,col)   = Rneg(k);
                    lag_map(row,col) = lag_vals(li) - p;
                end
            end
        end
    end

    if opts.showwaitbar && ishandle(hWB), waitbar(p/maxSteps, hWB); end
end
if opts.showwaitbar && exist('hWB','var') && ishandle(hWB), close(hWB); end

%% ───────── Visualise the evolving regressors ─────────
if opts.plot
    % Trim unused cells
    evolvingReg = evolvingReg(1:regCounter);

    % Build common time-axis in seconds
    t = (0:numel(evolvingReg{1})-1) .* opts.TR ./ opts.interp_factor;

    figure('Name','Evolving regressors','Color','w'); hold on
    cmap = jet(numel(evolvingReg));
    for ii = 1:numel(evolvingReg)
        plot(t, zscore(double(evolvingReg{ii})), 'Color', cmap(ii,:));
    end
    xlabel('Time (s)'); ylabel('z-scored amplitude');
    title('Average time-course of each updated regressor');
    colormap(cmap); cb = colorbar;
    cb.Label.String = 'Iteration (cool → warm)';
    box on, grid on
end

%% ---------------- Max‑r extraction ---------------------------
[maxR, bestCol] = max(r_map,[],2);
lag_best = NaN(xx*yy*zz,1,'single');
for k = 1:length(coordinates)
    lag_best(coordinates(k)) = lag_map(coordinates(k), bestCol(coordinates(k)));
end

lag_sec = opts.TR * (lag_best / opts.interp_factor);

%% ---------------- Output maps -------------------------------
lag_img = mask .* reshape(lag_sec, [xx yy zz]);
r_map   = mask .* reshape(maxR  , [xx yy zz]);

maps.recursiveLag_map = lag_img;
maps.recursiveR_map  = r_map;

if opts.smoothmap
    temp = opts.FWHM;
    opts.FWHM = [5 5 5];
    maps.recursiveLag_map = filterData(single(maps.recursiveLag_map), mask, mask, opts);
    opts.FWHM = temp; clear temp
end

%% ---------------- Save maps -------------------------------

savedir = opts.recursivedir
saveMap(cast(mask.*maps.recursiveLag_map,opts.mapDatatype), savedir, 'recursive_lag_map', opts.info.map, opts);
saveMap(cast(mask.*r_map,opts.mapDatatype), savedir, 'recursive_r_map', opts.info.map, opts);

end
