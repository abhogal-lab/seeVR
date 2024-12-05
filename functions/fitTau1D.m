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

    if isfield(opts, 'maxTau'); else; opts.maxTau = 300; end % Max exponential dispersion time constant

    %t = double(opts.TR:opts.TR:length(probe) * opts.TR);
    t = 1:1:length(probe);

    options = optimoptions('lsqcurvefit', 'Display', 'none', 'FunctionTolerance', 1.0000e-8, ...
        'StepTolerance', 1.0000e-8, 'MaxIter', 150);

    nr_params = 3;

    % Initialize weights: zero weight for zero values in ts, otherwise weight = 1
    weights = ts ~= 0;

    model = @(a, t) a(1) * rescale(real(ifft(ifftshift(fftinput(probe) .* fftexponential(a(2), t))))) + a(nr_params);
    a0 = [mean(ts(weights)) / mean(probe), 1, mean(ts(weights))];
    lb = [0, 0, -Inf];
    ub = [Inf, opts.maxTau, Inf];

    % Weighted least squares residual function
    weighted_residual = @(a) weights .* (ts - model(a, t));

    tic
    try
        % Use lsqnonlin for custom residual calculation
        b = lsqnonlin(weighted_residual, a0, lb, ub, options);
    catch
        % Continue on errors
        disp('Error: please check inputs');
    end

    % Generate fitted values
    y_fit = model(b, t);

    % Plot the results
    figure;
    plot(t, ts, 'r-', t, y_fit, 'b-');
    legend('Observed Data', 'Fitted Model');
    title('Observed vs. Fitted Model');
end



