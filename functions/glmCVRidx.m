% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <glmCVRidx: generates internally normalized CVR maps (i.e. without using respiratory data) >
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
% this function is inspited by the following manuscript:
% Cerebrovascular Reactivity Mapping Using Resting-State BOLD Functional MRI in Healthy Adults and Patients with Moyamoya Disease
% https://doi.org/10.1148/radiol.2021203568
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6494444/
% glmCVRidx uses a general linear model to evaluate the voxel-wise signal
% response relative to some reference signal. Specific data frequency bands
% can be isolated (i.e. for resting state fMRI analysis).
%
% data: input timeseries data (i.e. 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% refmask: binary mask that defines the reference signals of interest
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.fpass opts.headers.map, opts.info.map, opts,niiwrite, opts.resultsdir
%
% CVRidx_map: map of beta values (normalized to reference betas if
% opts.normalizeCVRidx = 1)
%
% BP_ref: the bandpassed reference signal (defined by opts.fpass) from the
% ROI defined by refmask
%
% bpData: a bandpassed version of the input data
function [CVRidx_map BP_ref bpData] = glmCVRidx(data, mask, refmask , opts)

warning('off');
global opts;

if isfield(opts,'fpass'); else; opts.fpass = [0.000001 0.08]; end  %default frequency band
if isfield(opts,'smoothmap'); else; opts.fpass = 0; end  %default is to not smooth output maps
if isfield(opts,'niiwrite'); else; opts.niiwrite = 0; end %depending on how data is loaded this can be set to 1 to use native load/save functions
if ispc
    opts.CVRidxdir = [opts.resultsdir,'CVRidx\']; mkdir(opts.CVRidxdir);
else
    opts.CVRidxdir = [opts.resultsdir,'CVRidx/']; mkdir(opts.CVRidxdir);
end

data = double(data);
[xx yy zz N] = size(data);


%get voxels
[voxels coordinates] = grabTimeseries(data, mask);
[ref_voxels ref_coordinates] = grabTimeseries(data, refmask);

Fs = 1/opts.TR;
t = 0:1/Fs:1-1/Fs;
N = size(voxels,2);
dF = Fs/N;
mean_ts = meanTimeseries(data,refmask);
freq = 0:Fs/length(mean_ts):Fs/2;
Lowf = opts.fpass(1); Highf = opts.fpass(2); %Hz

%check whether signal processing toolbox is available
if license('test','signal_toolbox') == 1
% if opts.filter_order is specified
if isfield(opts,'filter_order')
    [b,a] = butter(opts.filter_order,2*[Lowf, Highf]/Fs);
    BP_ref = filtfilt(b,a,(mean_ts));
else
    % else find optimal filter (automated)
    SSres = []; BP = [];
    for forder = 1:5
        [b,a] = butter(forder,2*[Lowf, Highf]/Fs);
        BP(forder,:) = filtfilt(b,a,(mean_ts));
        SSres(forder,:)= sum((mean_ts - BP(forder,:).^2),2);
    end
    %select closest filter
    [M,I] = min(SSres);
    [b,a] = butter(I,2*[Lowf, Highf]/Fs);
    BP_ref = filtfilt(b,a,(mean_ts));
    opts.filter_order = I;
end

%bandpass signals of interest contained within mask
BP_V = zeros(size(voxels));
parfor ii=1:size(BP_V,1)
    BP_V(ii,:) = filtfilt(b,a,voxels(ii,:));
end

%bandpass reference signals contained within refmask
BP_rV = zeros(size(ref_voxels));
parfor ii=1:size(BP_rV,1)
    BP_rV(ii,:) = filtfilt(b,a,ref_voxels(ii,:));
end

else
isplot = 0;
BP_ref = bpfilt(mean_ts, Lowf, Highf, isplot);
BP_ref = BP_ref';
%bandpass signals of interest contained within mask
BP_V = zeros(size(voxels));
BP_V = bpfilt(voxels,Lowf, Highf, isplot);

%bandpass reference signals contained within refmask
BP_rV = zeros(size(ref_voxels));
BP_rV= bpfilt(ref_voxels, Lowf, Highf, isplot);

end

figure; hold on
plot(rescale(mean_ts));
plot(rescale(BP_ref+mean(mean_ts)));
title('reference regressor before and after bandpass');
saveas(gcf,[opts.figdir,'reference_regressor.fig']);
xlabel('image volumes')
ylabel('a.u.')
legend('mean reference time-series', 'BP reference mean-timeseries')

%regress whole brain & reference signals against band passed reference
BP_ref = rescale(BP_ref); %reference regressor rescaled
D = [ones([length(BP_ref) 1]) BP_ref'];
coef = D\BP_V';
rcoef = D\BP_rV';
%reference average beta value
avg = mean(rcoef(2,:));

CVRidx = coef(2,:);
bpData = zeros(size(data));
bpData = reshape(bpData,[xx*yy*zz N]);
bpData(coordinates, :) = BP_V;
bpData = reshape(bpData,size(data));

CVRidx_map = zeros([1 numel(mask)]);
nCVRidx_map = zeros([1 numel(mask)]);

nCVRidx_map(1, coordinates) = CVRidx/avg; %normalized to reference response
CVRidx_map(1, coordinates) = CVRidx;

CVRidx_map = reshape(CVRidx_map, size(mask));
nCVRidx_map = reshape(nCVRidx_map, size(mask));

if opts.smoothmap
    % deprecated smoothing
    %CVRidx_map = smthData(CVRidx_map,mask,opts);
    %nCVRidx_map = smthData(nCVRidx_map,mask,opts);
    CVRidx_map = filterData(CVRidx_map,data(:,:,:,1),mask,opts);
    nCVRidx_map =filterData(nCVRidx_map,data(:,:,:,1),mask,opts);
end

nCVRidx_map(nCVRidx_map > 30) = 0;
nCVRidx_map(nCVRidx_map < -30) = 0;
CVRidx_map(CVRidx_map > 30) = 0;
CVRidx_map(CVRidx_map < -30) = 0;

if opts.niiwrite
    cd(opts.CVRidxdir);
    niftiwrite(nCVRidx_map,'normCVRidx_map',opts.info.map);
else
    saveImageData(nCVRidx_map,opts.headers.map,opts.CVRidxdir,'normCVRidx_map.nii.gz',64);
end

%remove strange values
nCVRidx_map(nCVRidx_map > 30) = 0;
nCVRidx_map(nCVRidx_map < -30) = 0;

if opts.niiwrite
    cd(opts.CVRidxdir);
    niftiwrite(CVRidx_map,'CVRidx_map',opts.info.map);
else
    saveImageData(CVRidx_map,opts.headers.map,opts.CVRidxdir,'CVRidx_map.nii.gz',64);
end
nCVRidx_map(nCVRidx_map == 0) = NaN;
CVRidx_map(CVRidx_map == 0) = NaN;
end

