% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <fitTau1D: calculates the dispersion time constant between in input
% signal and an associated MRI signal timeseries
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
%
% probe: usually PetCO2/O2 trace
%
% ts:    input timeseries data (i.e. average BOLD signal over an ROI)
%
% OPTIONS (opts)
%   opts.TR        : repetition time [s]  (required)
%   opts.maxTau    : max tau (default 300)
%   opts.maxDelay  : max |delay| in seconds (default 5)
%   opts.useDelay  : if true  -> fit delay (4 params: scale,tau,offset,delay)
%                    if false -> no delay term (3 params: scale,tau,offset)
%
function [y_fit, b] = fitTau1D(probe, ts, opts)

    global opts   %#ok<GVMIS>  % keep for backwards compatibility

    ts = double(ts);
    if iscolumn(probe); probe = probe'; end
    if iscolumn(ts);    ts    = ts';    end

    % ---- defaults ----
    if ~isfield(opts,'maxTau');    opts.maxTau    = 300; end
    if ~isfield(opts,'maxDelay');  opts.maxDelay  = 5;   end
    if ~isfield(opts,'useDelay');  opts.useDelay  = false; end

    if ~isfield(opts,'TR')
        error('fitTau1D:MissingTR','opts.TR must be provided (TR in seconds)');
    end

    % Time vector in seconds (as in your original)
    t = double(opts.TR : opts.TR : length(probe)*opts.TR);

    options = optimoptions('lsqcurvefit', ...
        'Display','none', ...
        'FunctionTolerance',1e-8, ...
        'StepTolerance',1e-8, ...
        'MaxIter',150);

    % === MODEL DEFINITIONS =================================================
    %
    % Base dispersion model (no delay), 3 params:
    %   a = [scale, tau, offset]
    baseModel = @(a, tLocal, probeLocal) ...
        a(1) * rescale(real(ifft(ifftshift( ...
               fftinput(probeLocal) .* fftexponential(a(2), a(1), tLocal))))) ...
        + a(3);

    % If delay is enabled, we will shift the probe (in samples) BEFORE
    % applying the same dispersion model.
    %
    % Full model with delay, 4 params:
    %   a = [scale, tau, offset, delay(sec)]
    delayedModel = @(a, tLocal, probeLocal, TR) ...
        baseModel(a(1:3), tLocal, shiftProbe(probeLocal, a(4), TR));
    % =======================================================================

    % === PARAMETER SETUP ===================================================
    if ~opts.useDelay
        % ---- 3-parameter fit: [scale, tau, offset] ----
        nr_params = 3;

        % Initial guess as in your original spirit
        a0 = [mean(ts)/mean(probe), 1, mean(ts)];
        lb = [0,   0,       -Inf];
        ub = [Inf, opts.maxTau, Inf];

        % You were actually using [range(ts) 20 0] as starting point;
        % to keep that behaviour:
        startPars = [range(ts), 20, 0];

        model = @(a,tLocal) baseModel(a, tLocal, probe);

    else
        % ---- 4-parameter fit: [scale, tau, offset, delay] ----
        nr_params = 4;

        a0 = [mean(ts)/mean(probe), 1, mean(ts), 0];
        lb = [0,   0,       -Inf,         -opts.maxDelay];
        ub = [Inf, opts.maxTau,  Inf,      opts.maxDelay];

        startPars = [range(ts), 20, 0, 0];

        model = @(a,tLocal) delayedModel(a, tLocal, probe, opts.TR);
    end

    % Dummy preallocation (not strictly needed but kept from original)
    b = nan(1, nr_params);

    % === FIT ===============================================================
    try
        b = lsqcurvefit(model, startPars, t, ts, lb, ub, options);
    catch
        disp('error; please check inputs');
        b = startPars;
    end

    % Evaluate final fit
    y_fit = model(b, t);

    % Extract scale & tau
    scaleFit = b(1);
    tauFit   = b(2);

    fprintf('fitTau1D: scale = %.4g, tau = %.4g\n', scaleFit, tauFit);

    % === BUILD SCALED PROBE FOR DISPLAY ====================================
    if opts.useDelay
        % Use the delayed probe actually used in the model
        probe_used = shiftProbe(probe, b(4), opts.TR);
    else
        probe_used = probe;
    end

    % Scale probe to similar amplitude range as ts/model
    minVal = min([ts(:); y_fit(:)]);
    maxVal = max([ts(:); y_fit(:)]);
    probe_scaled = rescale(probe_used, minVal, maxVal);

    % === PLOT ==============================================================
    figure;
    hold on;
    plot(t, ts,          'r-',  'LineWidth',1.4); % Observed data
    plot(t, y_fit,       'b-',  'LineWidth',1.4); % Fitted model
    plot(t, probe_scaled,'g--', 'LineWidth',1.2); % Scaled (and delayed) probe
    hold off;
    legend('Observed Data','Fitted Model','Scaled Probe','Location','best');
    xlabel('Time (s)');
    ylabel('Signal (a.u.)');

    if opts.useDelay
        title(sprintf('Observed vs. Fitted Model  (scale = %.3g, tau = %.3g, delay = %.3g s)', ...
              scaleFit, tauFit, b(4)));
    else
        title(sprintf('Observed vs. Fitted Model  (scale = %.3g, tau = %.3g, delay = 0)', ...
              scaleFit, tauFit));
    end
end

% === HELPER: shift probe by delay (seconds) using circular shift ==========
function pShift = shiftProbe(probe, delaySec, TR)
    % delaySec > 0  → probe is shifted forward in time (lags ts)
    % delaySec < 0  → probe leads
    shiftSamples = round(delaySec / TR);
    if shiftSamples == 0
        pShift = probe;
    else
        pShift = circshift(probe, shiftSamples);
    end
end
