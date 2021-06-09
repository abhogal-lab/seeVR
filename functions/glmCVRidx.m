function [CVRidx BP_ref bpBOLD] = glmCVRidx(data, mask, refmask , opts)
%Script written by ALEX BHOGAL a.bhogal@umcutrecht.nl
%resting state BOLD or scrubbed BOLD is used to generate resting state CVR
%maps based on band-passed signals
%data: timeseries data
%mask: mask of voxels of interest (normally whole-brain mask)
%refmask: mask of region from which regressor signals of interest are
%calculated. This can be whole-brain or cerebellum or sinus (for example)
%see: Cerebrovascular Reactivity Mapping Using Resting-State BOLD Functional MRI in Healthy Adults and Patients with Moyamoya Disease
% https://doi.org/10.1148/radiol.2021203568
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6494444/

global opts;
opts.CVRidxdir = [opts.resultsdir,'CVRidx\']; mkdir(opts.CVRidxdir);

data = double(data);
[xx yy zz N] = size(data);

opts.voxelsize = opts.headers.map.dime.pixdim(2:4)

%get voxels
[voxels coordinates] = grabTimeseries(data, mask);
[ref_voxels ref_coordinates] = grabTimeseries(data, refmask);
dim = opts.headers.map.dime.pixdim(2:4)
Fs = 1/opts.TR;
t = 0:1/Fs:1-1/Fs;
N = size(voxels,2);
dF = Fs/N;
mean_ts = meanTimeseries(data,refmask);
freq = 0:Fs/length(mean_ts):Fs/2;
Lowf = opts.fpass(1); Highf = opts.fpass(2); %Hz

% filter properties (automated)
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
filter_order = I

figure; plot(rescale(mean_ts)); hold on; plot(rescale(BP_ref+mean(mean_ts))); title('reference regressor');
saveas(gcf,[opts.figdir,'reference_regressor.fig']);

%bandpass whole brain signals
BP_V = zeros(size(voxels));
parfor ii=1:size(BP_V,1)
    BP_V(ii,:) = filtfilt(b,a,voxels(ii,:));
end
%bandpass reference signals
BP_rV = zeros(size(ref_voxels));
parfor ii=1:size(BP_rV,1)
    BP_rV(ii,:) = filtfilt(b,a,ref_voxels(ii,:));
end

%regress whole brain & reference signals against band passed reference
BP_ref = rescale(BP_ref); %reference regressor rescaled
D = [ones([length(BP_ref) 1]) BP_ref'];
coef = D\BP_V';
rcoef = D\BP_rV';
%reference average
avg = mean(rcoef(2,:));

CVRidx = coef(2,:);
bpBOLD = zeros(size(data));
bpBOLD = reshape(bpBOLD,[xx*yy*zz N]);
bpBOLD(coordinates, :) = BP_V;
bpBOLD = reshape(bpBOLD,size(data));

CVRidx_map = zeros([1 numel(mask)]);
if opts.normalizeCVRidx
    CVRidx_map(1, coordinates) = CVRidx/avg; %normalized to reference response
else
    CVRidx_map(1, coordinates) = CVRidx;
end
CVRidx_map = reshape(CVRidx_map, size(mask));

if opts.smoothmap
    CVRidx_map = smthData(CVRidx_map,mask,opts);
end
if opts.normalizeCVRidx
    saveImageData(CVRidx_map,opts.headers.map,opts.CVRidxdir,'normCVRidx_map.nii.gz',64)
else
    saveImageData(CVRidx_map,opts.headers.map,opts.CVRidxdir,'CVRidx_map.nii.gz',64)
end

end

