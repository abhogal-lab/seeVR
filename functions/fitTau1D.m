% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht
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
% ts: input timeseries data (i.e. average BOLD signal calculated over an ROI using meanTimeseries)
% 
% probe: usually PetCO2/O2 trace
%
% weighting: Added a weighting such that when input data contains volumes
% (or timeseries) with zero values, those values are weighted less in the
% fit. This means that rising or training edges can be isolated
% independently
%
function [y_fit, b] = fitTau1D(probe, ts, opts)
    global opts
    ts = double(ts);
    if iscolumn(probe); probe = probe'; end
    if iscolumn(ts); ts = ts'; end

    if isfield(opts, 'maxTau'); else; opts.maxTau = 300; end
    if isfield(opts, 'maxDelay'); else; opts.maxDelay = 5; end

    t = 1:1:length(probe); N = numel(t);
    fs = 1 / opts.TR;
    fvec = (0:N-1) * fs / N;

    options = optimoptions('lsqnonlin', 'Display', 'none', 'FunctionTolerance', 1e-6, ...
        'StepTolerance', 1e-6, 'MaxIter', 150);

    weights = double(ts ~= 0);

fvec = fftshift( (-floor(N/2):ceil(N/2)-1) * fs / N ); % Symmetric around zero
probe_fft = fftshift(fft(rescale(probe)));

phaseShift = @(d) exp(-1i * 2 * pi * fvec * d);

model = @(a, t) real(ifft(ifftshift(probe_fft .* phaseShift(a(4)) .* ...
                                    fftshift(fftexponential(a(2), a(1), t))))) + a(3);

    a0 = [mean(ts(weights)) / mean(probe), 1, mean(ts(weights)), 0];
    lb = [-Inf, 0, -Inf, -opts.maxDelay];
    ub = [Inf, opts.maxTau, Inf, opts.maxDelay];

    weighted_residual = @(a) weights .* (ts - model(a, t));

    try
        b = lsqnonlin(weighted_residual, a0, lb, ub, options);
    catch
        disp('Error: please check inputs');
        b = a0;
    end

    y_fit = model(b, t);

    figure;
    plot(t, ts, 'r-', t, y_fit, 'b-');
    legend('Observed Data', 'Fitted Model');
    title('Observed vs. Fitted Model');
end
