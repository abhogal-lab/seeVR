% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht
% a.bhogal@umcutrecht.nl
% <genGS: GLM based approach to regress out nuisance and data signals in provide residual timeseries data >
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

function [cleanData, res_ts, resData] = genGS(data, mask, nuisance, probe, opts)
% genGS — GLM-based residualization to make a "pseudo resting-state" scan
%
% Inputs:
%   data     : 3D/4D fMRI (x,y,z[,t])
%   mask     : binary ROI mask (x,y,z)
%   nuisance : T×K_n nuisance regressors (or []), T = #timepoints
%   probe    : T×K_p main probe(s) (or [])
%   opts     : options structure
%
% Main opts:
%   opts.residual_mode
%       'centered' : cleanData = centered residuals
%       'baseline' : cleanData = residuals + voxel-wise baseline
%
%   opts.baseline_method
%       'global_mean'   : use voxel-wise mean over full time series
%       'first_percent' : use voxel-wise mean over first X percent of time points
%       'first_n'       : use voxel-wise mean over first N time points
%
%   opts.baseline_percent
%       Used when opts.baseline_method = 'first_percent'
%       Default = 10
%
%   opts.baseline_n
%       Used when opts.baseline_method = 'first_n'
%       Default = round(0.1*T)
%
%   opts.save_output
%       0/default : do not save cleanData
%       1         : save cleanData
%
% Outputs:
%   cleanData : final cleaned pseudo-resting-state data, determined by
%               opts.residual_mode.
%
%               If opts.residual_mode = 'centered':
%                   cleanData = Y - Yhat
%
%               If opts.residual_mode = 'baseline':
%                   cleanData = Y - Yhat + voxel-wise baseline
%
%   res_ts    : mean cleanData time series within mask
%
%   resData   : zero-centered residual volume time series, always:
%                   resData = Y - Yhat
%
% Notes:
% - resData is retained as a secondary/debug/legacy output.
% - cleanData is the primary output to use and save.
% - The baseline is not taken from the GLM intercept because the intercept
%   can become numerically unstable if regressors are not centered, scaled,
%   or are poorly conditioned.

warning('off');  % quiet rank warnings from \ if any

% ---- Defaults -------------------------------------------------------------
if nargin < 5 || isempty(opts)
    opts = struct;
end

% New preferred option
if ~isfield(opts,'residual_mode') || isempty(opts.residual_mode)

    % Backward compatibility with older option
    if isfield(opts,'keep_offset') && opts.keep_offset
        opts.residual_mode = 'baseline';
    else
        opts.residual_mode = 'centered';
    end
end

% Backward compatibility with older alias
if isfield(opts,'add_offset_to_residual') && opts.add_offset_to_residual
    opts.residual_mode = 'baseline';
end

% New preferred save option
if ~isfield(opts,'save_output') || isempty(opts.save_output)

    % Backward compatibility with older save option
    if isfield(opts,'save_cleaned')
        opts.save_output = opts.save_cleaned;
    else
        opts.save_output = 0;
    end
end

if ~isfield(opts,'baseline_method') || isempty(opts.baseline_method)
    opts.baseline_method = 'global_mean';
end

if ~isfield(opts,'baseline_percent') || isempty(opts.baseline_percent)
    opts.baseline_percent = 10;
end

if ~isfield(opts,'save_residual')
    opts.save_residual = 0;
end

% ---- Validate residual mode ----------------------------------------------
opts.residual_mode = lower(opts.residual_mode);

switch opts.residual_mode
    case {'centered','zero_centered','zero-centred','zero_centered_residual'}
        opts.residual_mode = 'centered';

    case {'baseline','offset','with_baseline','baseline_restored','mean_restored'}
        opts.residual_mode = 'baseline';

    otherwise
        error('Unknown opts.residual_mode: %s. Use ''centered'' or ''baseline''.', opts.residual_mode);
end

% ---- Basic checks & shapes ------------------------------------------------
originalWas3D = false;

if ndims(data) == 4
    [nx, ny, nz, T] = size(data);
else
    originalWas3D = true;
    [nx, ny, T] = size(data);
    nz = 1;
    data = reshape(data, nx, ny, 1, T);
end

mask = logical(mask);

% Allow 2D mask when data is nx × ny × T
if nz == 1 && isequal(size(mask), [nx, ny])
    mask = reshape(mask, nx, ny, 1);
end

if ~isequal(size(mask), [nx, ny, nz])
    error('Mask size must match spatial size of data.');
end

if ~any(mask(:))
    error('Mask is empty.');
end

% Ensure probe/nuisance are T×K or empty
if ~isempty(probe)
    if size(probe,1) < size(probe,2)
        probe = probe';
    end
    if size(probe,1) ~= T
        error('probe must have T rows.');
    end
end

if ~isempty(nuisance)
    if size(nuisance,1) < size(nuisance,2)
        nuisance = nuisance';
    end
    if size(nuisance,1) ~= T
        error('nuisance must have T rows.');
    end
end

% ---- Build combined regressor matrix -------------------------------------
if ~isempty(probe) && ~isempty(nuisance)
    Xr = [nuisance, probe];        % T × K
elseif ~isempty(probe)
    Xr = probe;                    % T × Kp
elseif ~isempty(nuisance)
    Xr = nuisance;                 % T × Kn
else
    Xr = [];                       % Only intercept
end

% ---- Mean-center regressors ----------------------------------------------
% This makes the GLM numerically better behaved. With an intercept included,
% mean-centering the regressors does not change the model space.
if ~isempty(Xr)
    Xr = Xr - mean(Xr, 1, 'omitnan');
    Xr(~isfinite(Xr)) = 0;
end

% ---- Extract time series from mask ---------------------------------------
[V_ts, coords] = grabTimeseries(data, mask);  % V × T
Y = V_ts.';                                   % T × V

% Remove bad voxels if present
badVox = any(~isfinite(Y), 1);
if any(badVox)
    warning('%d voxels contain NaN/Inf and will be set to zero in residual output.', sum(badVox));
    Y(:,badVox) = 0;
end

% ---- Voxel-wise baseline --------------------------------------------------
% This is only used when opts.residual_mode = 'baseline'.
% For dynamic scans, use 'first_percent' or 'first_n' rather than global_mean.
switch lower(opts.baseline_method)

    case {'global_mean','mean','full_mean'}
        baseline_idx = 1:T;

    case {'first_percent','first_pct','initial_percent'}
        pct = opts.baseline_percent;

        if pct <= 0 || pct > 100
            error('opts.baseline_percent must be > 0 and <= 100.');
        end

        baseline_n = max(1, round(T * pct / 100));
        baseline_idx = 1:baseline_n;

    case {'first_n','initial_n'}
        if ~isfield(opts,'baseline_n') || isempty(opts.baseline_n)
            opts.baseline_n = max(1, round(0.1 * T));
        end

        baseline_n = round(opts.baseline_n);

        if baseline_n < 1
            error('opts.baseline_n must be at least 1.');
        end

        baseline_n = min(baseline_n, T);
        baseline_idx = 1:baseline_n;

    otherwise
        error('Unknown opts.baseline_method: %s', opts.baseline_method);
end

voxelBaseline = mean(Y(baseline_idx,:), 1, 'omitnan');   % 1 × V

% Safety fallback if any voxel baseline is bad
badBase = ~isfinite(voxelBaseline);
if any(badBase)
    voxelBaseline(badBase) = mean(Y(:,badBase), 1, 'omitnan');
end

% Final safety fallback if full mean also failed
badBase = ~isfinite(voxelBaseline);
if any(badBase)
    voxelBaseline(badBase) = 0;
end

% ---- Design matrix with intercept ----------------------------------------
if isempty(Xr)
    D = ones(T,1);                 % T × 1
else
    D = [ones(T,1), Xr];           % T × p
end

% ---- Solve voxel-wise GLM -------------------------------------------------
coef = D \ Y;                      % p × V

% ---- Fitted signals and residuals ----------------------------------------
Yhat = D * coef;                   % T × V, includes intercept
E = Y - Yhat;                      % T × V, centered residuals

% ---- Create primary cleaned output ---------------------------------------
switch opts.residual_mode

    case 'centered'
        Eout = E;                  % T × V

    case 'baseline'
        Eout = E + voxelBaseline;  % T × V
end

% ---- Pack back into volumes ----------------------------------------------

% cleanData = primary cleaned pseudo-resting-state data
cleanData = zeros(nx*ny*nz, T);
cleanData(coords, :) = Eout.';     % V × T
cleanData = reshape(cleanData, [nx, ny, nz, T]);

% resData = always centered residual data
resData = zeros(nx*ny*nz, T);
resData(coords, :) = E.';          % V × T
resData = reshape(resData, [nx, ny, nz, T]);

% ---- Mean cleaned residual time series within mask ------------------------
res_ts = mean(Eout, 2, 'omitnan').';            % 1 × T

% ---- Return same shape as input if original was nx × ny × T ---------------
if originalWas3D
    cleanData = reshape(cleanData, [nx, ny, T]);
    resData   = reshape(resData,   [nx, ny, T]);
end

% ---- Optional plotting ----------------------------------------------------
try
    if isfield(opts,'resultsdir') && ~isempty(opts.resultsdir)
        figdir = opts.resultsdir;
    elseif isfield(opts,'figdir') && ~isempty(opts.figdir)
        figdir = opts.figdir;
    else
        figdir = pwd;
    end

    if ~exist(figdir, 'dir')
        mkdir(figdir);
    end

    limits = [1, T];

    f = figure('visible','on');
    tiledlayout(f, 4, 1, "Padding","compact","TileSpacing","compact");

    % Original mean
    nexttile;
    plot(mean(V_ts,1,'omitnan'),'k');
    xlim(limits);
    title('Original mean signal ROI');
    ylabel('a.u.')

    % Probe(s)
    if ~isempty(probe)
        nexttile;
        hold on
        for k = 1:size(probe,2)
            plot(probe(:,k));
        end
        xlim(limits);
        title('Probe regressor(s)');
        ylabel('a.u.')
    else
        nexttile;
        axis off;
        text(0.02,0.5,'No probe provided');
    end

    % Nuisance
    if ~isempty(nuisance)
        nexttile;
        hold on
        for k = 1:size(nuisance,2)
            plot(nuisance(:,k));
        end
        xlim(limits);
        title('Nuisance regressor(s)');
        ylabel('a.u.')
    else
        nexttile;
        axis off;
        text(0.02,0.5,'No nuisance provided');
    end

    % Cleaned residual mean
    nexttile;
    plot(res_ts,'k');
    xlim(limits);

    switch opts.residual_mode
        case 'centered'
            ttl = 'Cleaned residual mean ROI centered';

        case 'baseline'
            switch lower(opts.baseline_method)
                case {'global_mean','mean','full_mean'}
                    ttl = 'Cleaned residual mean ROI with global voxel-wise baseline restored';
                case {'first_percent','first_pct','initial_percent'}
                    ttl = sprintf('Cleaned residual mean ROI with first %.1f%% baseline restored', opts.baseline_percent);
                case {'first_n','initial_n'}
                    ttl = sprintf('Cleaned residual mean ROI with first %d-point baseline restored', numel(baseline_idx));
                otherwise
                    ttl = 'Cleaned residual mean ROI with voxel-wise baseline restored';
            end
    end

    title(ttl);
    xlabel('time');
    ylabel('a.u.')

    saveas(f, fullfile(figdir,'residual_signal.fig'));

catch
    % plotting is best-effort; ignore failures in headless/batch
end

% ---- Optional saving ------------------------------------------------------
% Preferred behavior: save only the primary output cleanData.
if opts.save_output
    saveMap(cleanData, opts.resultsdir, 'cleanData', opts.info.ts, opts);
end

% Legacy/debug option: optionally save zero-centered residuals separately.
if opts.save_residual
    saveMap(resData, opts.resultsdir, 'residualData', opts.info.ts, opts);
end

end