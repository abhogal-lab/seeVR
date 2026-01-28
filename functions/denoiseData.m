% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <denoiseData: performs wavelet-based denoising or temporal smoothing >
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
% Depending on the available toolboxes, this function will perform either a
% wavelet based de-noising or a temporal smoothing of timeseries data
%
% data: input timeseries data (i.e. 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% opts: options structure containing required variables for this specific
% function; i.e. for wavelets: opts.wdlevel, opts.family, opts.figdir.
% for temporal smoothing: opts.method, opts.winsize (see MATLAB smoothdata
% function)
%
% denData: a denoised/smoothed version of the input data

function denData = denoiseData(data, mask, opts)
% denoiseData: performs wavelet-based denoising (if available) or temporal smoothing
%
% data: 4D timeseries (X x Y x Z x T)
% mask: 3D binary mask (X x Y x Z)
% opts fields (optional):
%   opts.wavelet  (default=1)   : try wavelet denoising
%   opts.wdlevel  (default=3)   : wavelet decomposition level
%   opts.family   (default='db4'): wavelet family
%   opts.method   (default='movmean'): smoothdata method (fallback / non-wavelet path)
%   opts.winsize  (default=round(0.05*T)): smoothing window size
%   opts.verbose  (default=0)   : plot mean timeseries before/after
%   opts.TR       (required if verbose=1): repetition time (s)
%   opts.figdir   (optional)    : folder to save .fig
%
% NOTE:
% - This version guarantees denData is always assigned.
% - Avoids global opts (input opts is authoritative).
% - Falls back cleanly if Wavelet Toolbox (or required funcs) are missing.

% -------------------- defaults / checks --------------------
if nargin < 3 || isempty(opts), opts = struct(); end

if ~isfield(opts,'wavelet'), opts.wavelet = 1; end
if ~isfield(opts,'wdlevel'), opts.wdlevel = 3; end
if ~isfield(opts,'family'),  opts.family  = 'db4'; end
if ~isfield(opts,'method'),  opts.method  = 'movmean'; end
if ~isfield(opts,'verbose'), opts.verbose = 0; end

[x, y, z, t] = size(data);

mask = logical(mask);
if ~isequal(size(mask), [x y z])
    error('denoiseData:MaskSizeMismatch', ...
        'mask must be size [%d %d %d] to match data spatial dims.', x, y, z);
end

tf   = class(data);
data = double(data);

% default winsize depends on T
if ~isfield(opts,'winsize') || isempty(opts.winsize)
    opts.winsize = max(1, round(0.05 * t)); % 5 percent window, min 1
end

% -------------------- choose method --------------------
didWavelet = false;

if opts.wavelet
    % Check toolbox license + existence of denoising function(s).
    hasWaveletLicense = false;
    try
        hasWaveletLicense = license('test','wavelet_toolbox');
    catch
        hasWaveletLicense = false;
    end

    % If your wavDenoise is YOUR function, exist(...) will be true if it's on path.
    % If you instead rely on MATLAB's Wavelet Toolbox denoisers, adapt accordingly.
    hasWavDenoiseFunc = (exist('wavDenoise','file') == 2);

    if hasWaveletLicense && hasWavDenoiseFunc
        try
            denData   = wavDenoise(data, mask, opts);
            didWavelet = true;
        catch ME
            warning('denoiseData:WaveletFailed', ...
                'Wavelet path failed (%s). Falling back to smoothing.', ME.message);
            didWavelet = false;
        end
    else
        % No error thrown here: we intentionally fall back.
        didWavelet = false;
        if ~hasWaveletLicense
            warning('denoiseData:NoWaveletLicense', ...
                'Wavelet Toolbox license unavailable. Falling back to smoothing.');
        elseif ~hasWavDenoiseFunc
            warning('denoiseData:MissingWavDenoise', ...
                'wavDenoise.m not found on path. Falling back to smoothing.');
        end
    end
end

% -------------------- fallback / smoothing --------------------
if ~didWavelet
    % Grab and smooth voxel time series
    [voxel_ts, coordinates] = grabTimeseries(data, mask);

    % smoothdata operates along dimension 2 for [Nvox x T]
    voxel_ts = smoothdata(voxel_ts, 2, opts.method, opts.winsize);

    denData = reshape(data, [x*y*z, t]);
    denData(coordinates, :) = voxel_ts;
    denData = reshape(denData, [x, y, z, t]);
end

% -------------------- optional QC plot --------------------
if opts.verbose
    if ~isfield(opts,'TR') || isempty(opts.TR)
        warning('denoiseData:MissingTR', ...
            'opts.TR not provided; using sample index on x-axis for verbose plot.');
        xdata = 1:t;
        xlab  = 'time (frames)';
    else
        xdata = opts.TR:opts.TR:(opts.TR*t);
        xlab  = 'time (s)';
    end

    figure; hold on;
    plot(xdata, meanTimeseries(data,   mask), 'k', 'LineWidth', 2);
    plot(xdata, meanTimeseries(denData,mask), 'c', 'LineWidth', 2);
    title('noise removal');
    ylabel('absolute signal'); xlabel(xlab);
    legend('original','denoised/smoothed');
    set(gcf, 'Units', 'pixels', 'Position', [200, 500, 600, 160]);
    xlim([xdata(1) xdata(end)]);
    hold off

    % Save if requested
    figName = 'denoiseData_QC.fig';
    if isfield(opts,'figdir') && ~isempty(opts.figdir) && exist(opts.figdir,'dir')
        saveas(gcf, fullfile(opts.figdir, figName));
    else
        saveas(gcf, fullfile(pwd, figName));
    end
end

denData = cast(denData, tf);

end
