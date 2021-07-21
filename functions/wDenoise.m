%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [ wdenData ssimval ]  = wDenoise(data, mask, opts)
%written by Alex Bhogal a.bhogal@umcutrecht.nl
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
% WB coordinates
[x,y,z,dyn] = size(data);
data(data == 0) = NaN;
wdenData = nan([x*y*z dyn]);
[voxel_ts, coordinates] = grabTimeseries(data, mask);

fd = zeros(size(voxel_ts));

%denoise signals
fd = wdenoise(voxel_ts',opts.level,'Wavelet',opts.family,'DenoisingMethod',opts.DenMeth,'ThresholdRule',opts.ThreshRule,'NoiseEstimate',opts.NoEst);
fd = fd';

wdenData (coordinates, :)  = fd;
wdenData  = reshape(wdenData ,[ x y z dyn]);

% compute the structural similarity index to evaluate denoising
% disp('computing structural similarity index')
% for ii=1:dyn
% ssimval(ii,:) = ssim(wdenBOLD(:,:,:,ii),data(:,:,:,ii));
% end
ssimval = [];
disp(['finished discreet wavelet denoising in: ',int2str(cputime-t),' seconds'])

end