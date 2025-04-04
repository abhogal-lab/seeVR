% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <bandpassfilt: performs a bandpass filter on input data >
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
function [bpData, bpMeanSignal] = bandpassfilt(data, mask, opts)
% bandpassfilt performs a bandpass filter on 4D input fMRI data within a given mask.
%
% Inputs:
% - data: 4D BOLD MRI dataset [X,Y,Z,T]
% - mask: binary ROI mask defining voxels to filter
% - opts: options structure with fields:
%       opts.TR - repetition time (s)
%       opts.fpass - frequency pass-band [low high] (Hz)
%       opts.filter_order (optional) - order of Butterworth filter
%
% Outputs:
% - bpData: bandpass-filtered dataset
% - bpMeanSignal: mean filtered signal within the mask
global opts
warning('off');
tic;

if ~isfield(opts,'fpass')
    opts.fpass = [0.000001 0.1164]; % default band: approx. 0.000001 - 0.1164 Hz
end

% Data dimensions
[xx, yy, zz, N] = size(data);
data = double(data);
mask = double(mask);
% Extract voxels within mask

[voxels_ts, coordinates] = grabTimeseries(data, mask);

% Mean ROI time series
mean_ts = mean(voxels_ts, 1);
Fs = 1 / opts.TR;
Lowf = opts.fpass(1);
Highf = opts.fpass(2);

% Define Butterworth filter
if license('test','signal_toolbox')
    if isfield(opts,'filter_order')
        [b,a] = butter(opts.filter_order, 2*[Lowf, Highf]/Fs);
    else
        % Automatically select optimal filter order (1 to 4)
        SSres = zeros(1,4);
        for order = 1:4
            [b_temp,a_temp] = butter(order, 2*[Lowf, Highf]/Fs);
            BP_temp = filtfilt(b_temp,a_temp,mean_ts);
            SSres(order) = sum((mean_ts - BP_temp).^2);
        end
        [~, best_order] = min(SSres);
        [b,a] = butter(best_order, 2*[Lowf, Highf]/Fs);
        opts.filter_order = best_order;
    end

    % Bandpass filter voxel time series
    BP_V = zeros(size(voxels_ts));
    parfor ii = 1:size(voxels_ts,1)
        try
        BP_V(ii,:) = filtfilt(b,a,voxels_ts(ii,:));
        catch
            disp(['error voxel :', int2str(ii)])
        BP_V(ii,:) = zeros(size(voxels_ts,2));
        end
    end
else
    error('Signal Processing Toolbox is required for Butterworth filtering.');
end

% Reconstruct filtered data
bpData = zeros(xx*yy*zz, N);
bpData(coordinates,:) = BP_V;
bpData = reshape(bpData, [xx, yy, zz, N]);

% Output mean filtered signal within mask
bpMeanSignal = filtfilt(b,a,mean(voxels_ts,1));

elapsedTime = toc;
disp(['Finished bandpass filtering in ', num2str(elapsedTime), ' seconds.']);

end