% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
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
% *************************************************************************
% This function uses input nuissance regressors and data regressors
% (i.e. response models) to remove input data of all explainable
% signal responses. The result is a residual 'global' signal that can be
% added as a further nuisance regressor (see scrubData) or used as 'pseudo
% resting-state' data.
%
% data: input timeseries data (i.e. 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% nuisance: and array of nuisance regressors (or a single regressor) having
% the same number of time-points as the input data.
%
% probe: an array of data-probes (explanatory variables) having the same
% number of time-points as the input data
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.figdir
%
% cleanData: residual data where nuisance and probe regressor explained
% signal variance has been removed. This data retains the linear offset
% from the GLM regression; in this way genGS can also be used to clean data
% of unwanted signal contributions (see scrubData)
%
% res_ts: the mean timeseries (cleanData) signal calculated within the ROI defined by
% the input mask
%
% resData: residual data after nuisance and probe regressor explained
% signal variance has been removed. Here the linear offset is also removed.

function [cleanData, res_ts, resData] = genGS(data, mask, nuisance, probe, opts)
% genGS2 — GLM-based residualization to make a "pseudo resting-state" scan
% Inputs:
%   data     : 3D/4D fMRI (x,y,z[,t])
%   mask     : binary ROI mask (x,y,z)
%   nuisance : T×K_n nuisance regressors (or []), T = #timepoints
%   probe    : T×K_p main probe(s) (or [])
%   opts     : struct (kept compatible with your original usage)
%
% Outputs:
%   cleanData : fitted data (intercept + modeled components) — same size as data
%   res_ts    : mean residual time series (no intercept) within mask
%   resData   : residual volume time series (no intercept) — same size as data
%
% Notes:
% - If both nuisance and probe are empty, this reduces to mean-removal (intercept model).
% - `cleanData` is the GLM fit; subtract it from `data` to get residuals+intercept.
%   We already return `resData` (residuals with intercept removed).

warning('off');  % quiet rank warnings from \ if any (ill-conditioning)
global opts;     % (requested: keep global opts aspect intact)
if ~isfield(opts,'save_cleaned'),   opts.save_cleaned = 0;   end

% ---- Basic checks & shapes ------------------------------------------------
if ndims(data) == 4
    [nx, ny, nz, T] = size(data);
else
    [nx, ny, T] = size(data); nz = 1; data = reshape(data, nx, ny, 1, T);
end
mask = logical(mask);
if any(size(mask) ~= [nx, ny, nz])
    error('Mask size must match spatial size of data.');
end

% Ensure probe/nuisance are T×K (or empty)
if ~isempty(probe)
    if size(probe,1) < size(probe,2), probe = probe'; end
    if size(probe,1) ~= T, error('probe must have T rows.'); end
end
if ~isempty(nuisance)
    if size(nuisance,1) < size(nuisance,2), nuisance = nuisance'; end
    if size(nuisance,1) ~= T, error('nuisance must have T rows.'); end
end

% Build combined regressor matrix (excluding intercept for now)
if ~isempty(probe) && ~isempty(nuisance)
    Xr = [nuisance, probe];        % T × K
elseif ~isempty(probe)
    Xr = probe;                    % T × Kp
elseif ~isempty(nuisance)
    Xr = nuisance;                 % T × Kn
else
    Xr = [];                       % Only intercept
end

% ---- Extract time series from mask ---------------------------------------
[V_ts, coords] = grabTimeseries(data, mask);  % V × T (V = #voxels in mask)
V = size(V_ts,1);
Y = V_ts.';                                    % T × V

% ---- Design matrix with intercept ----------------------------------------
if isempty(Xr)
    D = ones(T,1);                 % T × 1 (intercept only)
else
    D = [ones(T,1), Xr];          % T × p
end
p = size(D,2);

% ---- Solve voxel-wise GLM in one go --------------------------------------
% coef is p × V
coef = D \ Y;

% ---- Fitted signals and residuals ----------------------------------------
% Fitted (T×V):  Yhat = D * coef
Yhat = D * coef;

% Residuals (T×V): E = Y - Yhat
E = Y - Yhat;

% By definition, E has no intercept (the intercept is in Yhat).
% That makes E exactly what you want as "pseudo resting-state" residuals.

% ---- Pack back into volumes ----------------------------------------------
% Fitted -> cleanData
cleanData = zeros(nx*ny*nz, T);
cleanData(coords, :) = Yhat.';                     % V × T
cleanData = reshape(cleanData, [nx, ny, nz, T]);

% Residuals (no intercept) -> resData
resData = zeros(nx*ny*nz, T);
resData(coords, :) = E.';                          % V × T
resData = reshape(resData, [nx, ny, nz, T]);

% If original was 3D (single time series stack), return same shape
if nz == 1
    cleanData = reshape(cleanData, [nx, ny, T]);
    resData   = reshape(resData,   [nx, ny, T]);
end

% ---- Mean residual time series within mask --------------------------------
res_ts = meanTimeseries(resData, mask);   % 1 × T

% ---- Optional plotting (kept lightweight & robust) ------------------------
try
    if isfield(opts,'figdir') && ~isempty(opts.figdir)
        figdir = opts.figdir;
    else
        figdir = pwd;
    end
    limits = [1, T];

    f = figure('visible','on'); 
    tiledlayout(f, 4, 1, "Padding","compact","TileSpacing","compact");

    % Original mean
    nexttile; plot(meanTimeseries(data,mask),'k'); xlim(limits);
    title('Original mean signal (ROI)'); ylabel('a.u.')

    % Probe(s)
    if ~isempty(probe)
        nexttile; hold on
        for k = 1:size(probe,2), plot(probe(:,k)); end
        xlim(limits); title('Probe regressor(s)'); ylabel('a.u.')
    else
        nexttile; axis off; text(0.02,0.5,'No probe provided');
    end

    % Nuisance
    if ~isempty(nuisance)
        nexttile; hold on
        for k = 1:size(nuisance,2), plot(nuisance(:,k)); end
        xlim(limits); title('Nuisance regressor(s)'); ylabel('a.u.')
    else
        nexttile; axis off; text(0.02,0.5,'No nuisance provided');
    end

    % Residual mean
    nexttile; plot(res_ts,'k'); xlim(limits);
    title('Residual mean (ROI)'); xlabel('time'); ylabel('a.u.')

    saveas(f, fullfile(figdir,'residual_signal.fig'));

catch
    % plotting is best-effort; ignore failures in headless/batch
end
if opts.save_residual
    saveMap(resData, opts.resultsdir,'residualData',opts.info.ts, opts);
end

end
