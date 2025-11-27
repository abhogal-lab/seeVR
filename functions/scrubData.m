% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <scrubData: GLM based nuisance signal regression >
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
% Linear regression of input probes and nuisance regressors is performed on
% data in regions specified by the mask. Nuisance timeseries defined by
% nuisance regressors and co-efficients are summed and removed. For a
% similar function to determine residual signals after regression, see
% genGS.
%
% data: input timeseries data (i.e. 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% nuisance: and array of nuisance regressors (or a single regressor) having
% the same number of time-points as the input data. Generally these will be
% the rotations and translations derived from motion correction
%
% probe: an array of data-probes (explanatory variables) having the same
% number of time-points as the input data
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.figdir, opts.headers.ts, opts.resultsdir
%
% OUTPUTS
%
% cleanData: data scrubbed of nuisance signals
%
% np0: nuisance signals used during scrubbing
%
% np0: nuisance signals with high correlation to input that are rejected
%
function [cleanData, np0, np1] = scrubData2(data, mask, nuisance, probe, opts)
% scrubData2 — GLM-based nuisance removal (preserve probe effects)
%
% Inputs
%   data     : 3D/4D fMRI (x,y,z[,t])
%   mask     : logical ROI mask (x,y,z)
%   nuisance : T×K_n nuisance regressors (or []), T=#timepoints
%   probe    : T×K_p probe regressors used to *preserve* their variance (or [])
%   opts     : struct; fields used:
%              - niiwrite (0/1), save_cleaned (0/1), prep_nuisance (0/1),
%                disperse_probe (0/1), resultsdir, headers.ts / info.ts,
%                headers.map / info.map, figdir
%
% Outputs
%   cleanData : data with ONLY nuisance contribution removed
%   np0       : final nuisance regressors used
%   np1       : rejected/high-correlation nuisance regressors (from prepNuisance)

warning('off');
global opts;

% ---- Coerce types / shapes ----------------------------------------------
tf = class(data);
data = double(data);
mask = logical(mask);

if ndims(data) == 4
    [nx, ny, nz, T] = size(data);
else
    [nx, ny, T] = size(data); nz = 1;
    data = reshape(data, nx, ny, 1, T);
end
if any(size(mask) ~= [nx, ny, nz])
    error('Mask size must match data spatial size.');
end

% Ensure regressors are T×K (or empty)
if ~isempty(probe)
    probe = double(probe);
    if size(probe,1) < size(probe,2), probe = probe'; end
    if size(probe,1) ~= T, error('probe must have T rows.'); end
else
    probe = [];
end
if ~isempty(nuisance)
    nuisance = double(nuisance);
    if size(nuisance,1) < size(nuisance,2), nuisance = nuisance'; end
    if size(nuisance,1) ~= T, error('nuisance must have T rows.'); end
else
    nuisance = [];
end

% ---- Defaults ------------------------------------------------------------
if ~isfield(opts,'niiwrite'),       opts.niiwrite = 0;       end
if ~isfield(opts,'save_cleaned'),   opts.save_cleaned = 0;   end
if ~isfield(opts,'disperse_probe'), opts.disperse_probe = 0; end
if ~isfield(opts,'prep_nuisance'),  opts.prep_nuisance = 1;  end
if ~isfield(opts,'resultsdir'),     opts.resultsdir = pwd;   end

% ---- Prepare nuisance (filtering, orthogonalization, etc.) --------------
if opts.prep_nuisance && ~isempty(nuisance)
    disp('Preparing nuisance signals...');
    [np0, np1] = prepNuisance(nuisance, probe, opts);   % user-provided helper
else
    np0 = nuisance;
    np1 = [];
end

% ---- Optional probe dispersion (kept compatible with your code) ---------
if opts.disperse_probe && ~isempty(probe)
    disp('Adding dispersion terms to input probe');
    [~,~,probe] = convHRF(probe, opts);  % user-provided helper
    probe = probe';
end

% ---- Extract masked time series -----------------------------------------
[V_ts, coords] = grabTimeseries(data, mask);   % V × T
V = size(V_ts,1);
Y = V_ts.';                                    % T × V

% ---- Build design and remove ONLY nuisance contribution ------------------
% We use a full design with intercept + (optional) probe + nuisance,
% but we subtract only the fitted nuisance part: (nuisance * beta_nuis).
% This preserves probe effects and the intercept.
if isempty(np0) && isempty(probe)
    % Nothing to remove: return original data
    Y_clean = Y;
    beta_nuis = [];
elseif isempty(np0)
    % No nuisance -> nothing to remove; keep probe effects
    D = [ones(T,1), probe];            % T×(1+Kp)
    coef = D \ Y;                      % (1+Kp) × V
    Y_clean = Y;                       % do not remove probe
    beta_nuis = [];
else
    % With nuisance (and maybe probe)
    if isempty(probe)
        D = [ones(T,1), np0];          % T×(1+Kn)
        coef = D \ Y;                   % (1+Kn) × V
        beta_nuis = coef(2:end, :);     % Kn × V
        nuis_contrib = np0 * beta_nuis; % T×V
        Y_clean = Y - nuis_contrib;     % remove ONLY nuisance
    else
        D = [ones(T,1), probe, np0];    % T×(1+Kp+Kn)
        coef = D \ Y;                   % (1+Kp+Kn) × V
        % Nuisance betas are the LAST Kn rows (conditional on probe)
        Kn = size(np0,2);
        beta_nuis = coef(end-Kn+1:end, :);   % Kn × V
        nuis_contrib = np0 * beta_nuis;      % T×V
        Y_clean = Y - nuis_contrib;          % remove ONLY nuisance
    end
end

% ---- Pack cleaned data back to volume -----------------------------------
cleanData = zeros(nx*ny*nz, T);
cleanData(coords, :) = Y_clean.';               % V×T
cleanData = reshape(cleanData, [nx, ny, nz, T]);
cleanData = cast(cleanData, tf);

% ---- QC plots (best-effort) ---------------------------------------------
try
    limits = [1, T];
    if isfield(opts,'figdir') && ~isempty(opts.figdir)
        figdir = opts.figdir;
    else
        figdir = opts.resultsdir;
    end
    f = figure('visible','on');
    tiledlayout(f, 4, 1, "Padding","compact", "TileSpacing","compact");

    % Original mean
    nexttile; 
    plot(meanTimeseries(data,mask),'k'); 
    xlim(limits);
    title('Original mean (ROI)'); ylabel('a.u.')

    % Probe(s)
    nexttile;
    if ~isempty(probe)
        hold on; 
        for k=1:size(probe,2), plot(probe(:,k)); end
        hold off
        title('Probe regressor(s)'); xlim(limits);
    else
        axis off; text(0.02,0.5,'No probe provided');
    end

    % Nuisance
    nexttile;
    if ~isempty(np0)
        hold on; 
        for k=1:size(np0,2), plot(np0(:,k)); end
        hold off
        title('Nuisance regressor(s)'); xlim(limits);
    else
        axis off; text(0.02,0.5,'No nuisance provided');
    end

    % Cleaned mean (plot together with original)
    nexttile;
    plot(meanTimeseries(data,mask),'k','DisplayName','Original'); hold on
    plot(meanTimeseries(cleanData,mask),'r','DisplayName','Cleaned');
    xlim(limits);
    title('Original vs Cleaned mean (ROI)'); xlabel('time'); ylabel('a.u.')
    legend('show','Location','best'); hold off

    saveas(f, fullfile(figdir,'scrubData2.fig')); 
catch
    % ignore plotting errors in headless mode
end

% ---- Optional saving of cleaned data ------------------------------------
if opts.save_cleaned
    saveMap(cleanData, opts.resultsdir,'cleanData',opts.info.ts, opts);
end

% ---- QC metrics: Euclidean distance & MAPE (per voxel) -------------------
% Original masked series (double)
origY = V_ts.';            % T×V
cleanY= Y_clean;           % T×V
diffY = origY - cleanY;    % T×V (this is the removed nuisance component)

% Euclidean norm per voxel
euclid = sqrt(sum(diffY.^2, 1));          % 1×V

% MAPE per voxel (percent); avoid divide-by-zero using denom = max(|orig|, eps)
denom  = max(abs(origY), eps);
mape_v = mean(abs(diffY) ./ denom, 1) * 100;   % 1×V

% Map to volume
euclid_map = zeros(1, nx*ny*nz);
mape_map   = zeros(1, nx*ny*nz);
euclid_map(1, coords) = euclid;
mape_map(1,   coords) = mape_v;
euclid_map = reshape(euclid_map, [nx, ny, nz]);
mape_map   = reshape(mape_map,   [nx, ny, nz]);

% Save QC maps
saveMap(cast(euclid_map, opts.mapDatatype), opts.resultsdir,'euclidean_distance',opts.info.map, opts);
saveMap(cast(mape_map,   opts.mapDatatype), opts.resultsdir,'MAPE',opts.info.map, opts);

end
