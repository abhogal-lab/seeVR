%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [cleanData] = scrubData(data,mask, nuisance, probe, opts)
%written by Alex Bhogal a.bhogal@umcutrecht.nl
%Linear regression of input probes and nuisance regressors is performed on
%data in regions specified by the mask. Nuisance timeseries defined by
%nuisance regressors and co-efficients are summed and removed. For a
%similar function to determine residual signals after regression, see
%genGS.
warning('off');
global opts;

if isempty(probe); else
test1 = probe(1,:); test2 = probe(:,1);
if length(test1) > length(test2); probe = probe'; end; clear test1 test2
limits = [0 size(probe,1)];
probemap = plasma(size(probe,2)); %probemap = flip(probemap,1);
end
if isempty(nuisance); else
test1 = nuisance(1,:); test2 = nuisance(:,1);
if length(test1) > length(test2); nuisance = nuisance'; end; clear test1 test2
limits = [0 size(nuisance,1)];
nuisancemap = parula(size(nuisance,2)); nuisancemap = flip(nuisancemap,1);
end

[voxel_ts, coordinates] = grabTimeseries(data, mask);
[probe,~]=licols(probe);
probesize = size(probe,2);
np = [probe nuisance];


%perform regression
D = [ones([length(np) 1]) np];
np_coef = D\voxel_ts';

%nnintcp = squeeze(np_coef(1,:));
nuis = squeeze(np_coef(2+probesize:end,:));
nuis_TS = zeros(size(voxel_ts))';

parfor ii=1:size(nuisance,2)
    tmp = nuisance(:,ii)*nuis(ii,:);
    nuis_TS = nuis_TS+tmp;
end
switch ndims(data)
    case 4
        [x,y,z,dyn] = size(data);
        clean_data = voxel_ts - nuis_TS';
        cleanData = zeros([x*y*z,size(clean_data,2)]);
        cleanData(coordinates,:) = clean_data ;
        cleanData = reshape(cleanData,[x y z dyn]);
    case 3
        [x,y,dyn] = size(data);
        clean_data = voxel_ts - nuis_TS';
        cleanData = zeros([x*y,size(clean_data,2)]);
        cleanData(coordinates,:) = clean_data ;
        cleanData = reshape(cleanData,[x y dyn]);
end

figure;

if size(probe,2) > 1
subplot(5,1,1); hold on; for ii=1:size(probe,2); plot(np(:,ii), 'Color', probemap(ii,:)); end; title('data probe(s)'); xlim(limits);
else
subplot(5,1,1); hold on; for ii=1:size(probe,2); plot(np(:,ii), 'Color', 'k'); end; title('data probe(s)'); xlim(limits);
end
subplot(5,1,2); hold on; for ii=1:size(nuisance,2); plot(rescale(nuisance(:,ii)), 'Color', nuisancemap(ii,:)); end; title('rescaled nuisance regressors'); xlim(limits);
subplot(5,1,3); plot(meanTimeseries(data,mask),'k'); title('Original Data'); xlim(limits)
subplot(5,1,4); plot(nanmean(nuis_TS,2),'k'); title('nuisance Mean'); xlim(limits);
subplot(5,1,5); plot(meanTimeseries(cleanData,mask),'k'); title('Clean Data'); xlim(limits);

if isfield(opts,'figdir')
    saveas(gcf,[opts.figdir,'scrubData.fig']);
else
    saveas(gcf,[pwd,'\','scrubData.fig']);
end

if isfield(opts,'resultsdir')
else
    opts.resultsdir = [pwd,'\'];
end

if isfield(opts,'save_cleaned'); else opts.save_cleaned = 0; end
if opts.save_cleaned
saveImageData(cleanBOLD, opts.headers.ts, opts.resultsdir, 'cleanBOLD.nii.gz', 64); 
end


end

