function[ fSpec, vmeanSpec, avgSpec, freq ] = freqSpec(data, mask, opts)
% freqSpec calculates the frequency spectrum of a 4D timeseries for 
% voxels defined by a 3D mask.
%
% Outputs:
%   fSpec    - frequency spectrum per voxel (4D: x,y,z,freq)
%   vmeanSpec- average (mean) of voxel-wise spectra (i.e. mean over voxels)
%   avgSpec  - frequency spectrum of the average time series within the mask
%   freq     - vector of frequency values
%
% Inputs:
%   data     - 4D MRI timeseries data (x, y, z, time)
%   mask     - binary mask defining the region of interest (3D)
%   opts     - structure with fields:
%                TR          - repetition time
%                powerspec   - if 1, return power spectrum, else amplitude
%                showplots   - if 1, display plots
%                logtransform- if 1, apply log transformation to the plot
%                vLines      - (optional) vector of frequencies for vertical lines

% Set default options if not provided
global opts
if ~isfield(opts, 'powerspec');    opts.powerspec = 0; end
if ~isfield(opts, 'showplots');      opts.showplots = 0; end
if ~isfield(opts, 'logtransform');   opts.logtransform = 0; end

% Sampling frequency and Nyquist frequency
Fs = 1 / opts.TR;
max_freq = Fs/2;

% Set default vertical lines and filter out values beyond the Nyquist frequency
if ~isfield(opts, 'vLines')
    default_vLines = [0.08, 0.15, 0.2, 0.33];
    opts.vLines = default_vLines(default_vLines <= max_freq);
else
    % Ensure that provided vertical lines are within the valid frequency range
    opts.vLines = opts.vLines(opts.vLines <= max_freq);
end

% Extract the time series from voxels within the mask
[voxels, coordinates] = grabTimeseries(data, mask);
[x, y, z, N] = size(data);

% Frequency vector for the half-spectrum
freq = 0 : Fs/N : max_freq;

% --- Compute voxel-wise FFT ---
% FFT along the time dimension for each voxel
r_xdft = fft(voxels, [], 2);
% Keep only the positive half of the spectrum
r_xdft = r_xdft(:, 1:N/2+1);

% Convert FFT to either amplitude or power spectrum
if opts.powerspec
    r_psdx = (1/(Fs*N)) * abs(r_xdft).^2;
else
    r_psdx = (1/(Fs*N)) * abs(r_xdft);
end
% Adjust the spectrum except for DC and Nyquist components
r_psdx(2:end-1) = 2 * r_psdx(2:end-1);

% --- Compute average spectrum across voxels ---
% This is the mean over voxel-wise spectra
vmeanSpec = mean(r_psdx, 1);

% --- Compute spectrum of the average signal ---
% First, compute the average time series across voxels in the mask
averageTS = mean(voxels, 1);
% FFT of the average time series
avg_fft = fft(averageTS);
% Retain only the positive half of the spectrum
avg_fft = avg_fft(1:N/2+1);
if opts.powerspec
    avgSpec = (1/(Fs*N)) * abs(avg_fft).^2;
else
    avgSpec = (1/(Fs*N)) * abs(avg_fft);
end
avgSpec(2:end-1) = 2 * avgSpec(2:end-1);

% --- Plotting ---
if opts.showplots
    figure; set(gcf, 'color', 'w');
    hold on;
    
    % Plot the two main spectra curves
    if opts.logtransform
        h1 = plot(freq, log(vmeanSpec), 'LineWidth', 2, 'DisplayName', 'Mean voxel-wise spectrum (log)');
        %h2 = plot(freq, log(avgSpec), 'LineWidth', 2, 'DisplayName', 'Spectrum of average signal (log)');
    else
        h1 = plot(freq, vmeanSpec, 'LineWidth', 2, 'DisplayName', 'Mean voxel-wise spectrum');
       % h2 = plot(freq, avgSpec, 'LineWidth', 2, 'DisplayName', 'Spectrum of average signal');
    end

    % Generate distinct colors for each vertical line using the "lines" colormap
    nLines = numel(opts.vLines);
    lineColors = lines(nLines);
    
    % Plot vertical lines with custom colors and add them to the legend
    for i = 1:nLines
        xline(opts.vLines(i), 'LineStyle', '--', 'LineWidth', 1.5, 'Color', lineColors(i,:), ...
            'DisplayName', sprintf('%.3f Hz', opts.vLines(i)));
    end

    % Set title based on spectrum type
    if opts.powerspec
        if opts.logtransform
            title('Log of Power Spectrum');
        else
            title('Power Spectrum');
        end
    else
        if opts.logtransform
            title('Log of Amplitude Spectrum');
        else
            title('Amplitude Spectrum');
        end
    end

    legend('show');
    hold off;
end

% --- Assign per-voxel spectrum into a 4D matrix ---
% Create an empty array and fill in the voxels corresponding to the mask
fSpec = zeros([x*y*z, size(r_psdx,2)]);
fSpec(coordinates, :) = r_psdx;
fSpec = reshape(fSpec, [x, y, z, size(r_psdx,2)]);

disp('Returning spectra for the mask region (NOT log-transformed)');
saveas(gcf,fullfile(opts.figdir,'spectrum.png'))
saveas(gcf,fullfile(opts.figdir,'pectrum.fig'))
saveas(gcf,fullfile(opts.figdir,'spectrum.svg'))
end
