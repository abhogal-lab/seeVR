% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <wDenoise: a wrapper function to perform wavelet-based denoising using the MATLAB wavelet toolbox >
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
% This is a wrapper function which uses the wavelet toolbox to de-noise
% input signals. For more information see the MATLAB wdenoise function.
function [wdenData]  = wavDenoise(data, mask, opts)

global opts;
t = cputime;

%discrete wavelet transform ideal for denoising
%1. multi-level wavelet decomposition to obtain co-efficients (fine-scale
%analysys in sub-bands determined with HP and LP filters)
%2. identify suitable thresholding technique based on data
%3. threshold co-efficients and reconstruct the signal
if isfield(opts,'family'); else opts.family = 'db4'; end
if isfield(opts,'level'); else opts.level = 2; end
if isfield(opts,'DenMeth'); else opts.DenMeth = 'UniversalThreshold'; end
if isfield(opts,'ThreshRule'); else opts.ThreshRule = 'Hard'; end
if isfield(opts,'NoEst'); else opts.NoEst = 'LevelIndependent'; end
if isfield(opts,'verbose'); else; opts.verbose = 0; end %turn on/off select command output

% WB coordinates
[x,y,z,dyn] = size(data);
data(data == 0) = NaN;
wdenData = nan([x*y*z dyn]);
[voxel_ts, coordinates] = grabTimeseries(data, mask);

fd = zeros(size(voxel_ts));

%denoise signals
fd = wdenoise(double(voxel_ts'),opts.level,'Wavelet',opts.family,'DenoisingMethod',opts.DenMeth,'ThresholdRule',opts.ThreshRule,'NoiseEstimate',opts.NoEst);
fd = fd';

wdenData (coordinates, :)  = fd;
wdenData  = reshape(wdenData ,[ x y z dyn]);

disp(['finished discreet wavelet denoising in: ',int2str(cputime-t),' seconds']);

end