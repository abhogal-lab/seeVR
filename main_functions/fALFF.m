% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <fALFF: calculates (fractional) amplitude of low frequency fluctuations ((f)ALFF) >
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
function [ALFF_map, fALFF_map, rALFF_map, zALFF_map, zfALFF_map, zrALFF_map] = fALFF(data, mask, refmask, opts)
% fALFF - Computes ALFF, fALFF, and rALFF maps from 4D fMRI data.
%
% Usage:
%   [ALFF_map, fALFF_map, rALFF_map, zALFF_map, zfALFF_map, zrALFF_map] = ...
%       fALFF(data, mask, refmask, opts)
%
% Inputs:
%   data      - 4D fMRI data matrix (X × Y × Z × T), double or single precision.
%   mask      - 3D logical or binary mask defining voxels of interest.
%   refmask   - 3D logical or binary mask for the reference region (e.g., whole brain).
%   opts      - Struct containing processing options with required/optional fields:
%       .TR           - Repetition time in seconds (required).
%       .fpass        - 1x2 vector specifying frequency band [low high] in Hz (default: [0.01 0.15]).
%       .resultsdir   - Directory to save output maps (required).
%       .info         - Template NIfTI header info (required if opts.niiwrite = 1).
%       .headers.map  - Map header info (required if opts.niiwrite = 0).
%       .niiwrite     - Boolean flag to save as NIfTI (.nii) or custom (default: 0 = custom).
%
% Outputs:
%   ALFF_map   - Amplitude of low frequency fluctuations (raw ALFF).
%   fALFF_map  - ALFF normalized by full frequency spectrum (fractional ALFF).
%   rALFF_map  - ALFF normalized by reference region ALFF (relative ALFF).
%   zALFF_map  - Z-transformed ALFF map.
%   zfALFF_map - Z-transformed fALFF map.
%   zrALFF_map - Z-transformed rALFF map.
%
% Description:
%   The function computes three variants of ALFF:
%     1. ALFF: Raw amplitude in a low-frequency band (default 0.01–0.15 Hz).
%     2. fALFF: Fraction of ALFF relative to the full spectral amplitude.
%     3. rALFF: ALFF normalized to the mean ALFF of a reference region.
%
%   It also outputs z-scored versions of all three metrics. Optionally saves maps
%   as NIfTI or using a custom saving function (`saveImageData`).
%
% Example:
%   opts.TR = 2;
%   opts.resultsdir = './results';
%   opts.headers.map = header;
%   opts.niiwrite = 0;
%   [ALFF, fALFF, rALFF, zA, zf, zr] = fALFF(data, brainMask, brainMask, opts);
%
% References:
%   - Zang et al., 2007. Altered baseline brain activity in children with ADHD revealed by resting-state fMRI. Brain & Development.
%   - Zou et al., 2008. An improved approach to detection of amplitude of low-frequency fluctuation (ALFF) for resting-state fMRI. Front Syst Neurosci.
%   - Yang et al., 2007. Amplitude of low frequency fluctuation within visual areas revealed by resting-state functional MRI. NeuroImage.
%
% See also: fft, meanTimeseries, grabTimeseries, saveImageData, niftiwrite

global opts
data = double(data);
temp_info = opts.info.map;
temp_info.Datatype = 'double';

if isfield(opts,'niiwrite') == 0, opts.niiwrite = 0; end
if isfield(opts,'fpass') == 0, opts.fpass = [0.01 0.15]; end
opts.ALFFdir = fullfile(opts.resultsdir,'ALFF'); mkdir(opts.ALFFdir);

[xx, yy, zz, N] = size(data);
Fs = 1 / opts.TR;
dF = Fs / N;
freq = 0:Fs/N:Fs/2;
[Lowval, Lowidx] = min(abs(freq - opts.fpass(1)));
[Highval, Highidx] = min(abs(freq - opts.fpass(2)));

%% Compute reference region ALFF
refdata = meanTimeseries(data, refmask);
xdft = fft(refdata);
xdft = xdft(1:N/2+1);
psdx = (abs(xdft).^2)/N;
psdx(2:end-1) = 2*psdx(2:end-1);
refALFF = sum(sqrt(psdx(Lowidx:Highidx)));  % reference ALFF for rALFF

%% Voxel-wise ALFF
[voxels, coordinates] = grabTimeseries(data, mask);
xdft = fft(voxels, [], 2);
xdft = xdft(:, 1:N/2+1);
psdx = (abs(xdft).^2)/N;
psdx(:, 2:end-1) = 2 * psdx(:, 2:end-1);

vALFF = sum(sqrt(psdx(:, Lowidx:Highidx)), 2);  % raw ALFF
fvALFF = sum(sqrt(psdx), 2);                   % full spectrum

% Z-transformed values
zALFF = (vALFF - mean(vALFF)) / std(vALFF);
zfALFF = (vALFF ./ fvALFF - mean(vALFF ./ fvALFF)) / std(vALFF ./ fvALFF);
zrALFF = (vALFF / refALFF - mean(vALFF / refALFF)) / std(vALFF / refALFF);

%% Initialize output maps
mapSize = size(mask);
ALFF_map = zeros(mapSize);      zALFF_map = zeros(mapSize);
fALFF_map = zeros(mapSize);     zfALFF_map = zeros(mapSize);
rALFF_map = zeros(mapSize);     zrALFF_map = zeros(mapSize);

ALFF_map(coordinates)   = vALFF;
zALFF_map(coordinates)  = zALFF;

fALFF_map(coordinates)  = vALFF ./ fvALFF;
zfALFF_map(coordinates) = zfALFF;

rALFF_map(coordinates)  = vALFF / refALFF;
zrALFF_map(coordinates) = zrALFF;

% Reshape to 3D
ALFF_map   = reshape(ALFF_map, mapSize);
zALFF_map  = reshape(zALFF_map, mapSize);
fALFF_map  = reshape(fALFF_map, mapSize);
zfALFF_map = reshape(zfALFF_map, mapSize);
rALFF_map  = reshape(rALFF_map, mapSize);
zrALFF_map = reshape(zrALFF_map, mapSize);

%% Save results
saveMap(ALFF_map, 'ALFF_map', temp_info, opts);
saveMap(zALFF_map, 'zALFF_map', temp_info, opts);
saveMap(fALFF_map, 'fALFF_map', temp_info, opts);
saveMap(zfALFF_map, 'zfALFF_map', temp_info, opts);
saveMap(rALFF_map, 'rALFF_map', temp_info, opts);
saveMap(zrALFF_map, 'zrALFF_map', temp_info, opts);

end


