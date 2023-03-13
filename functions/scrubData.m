% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <scrubData: GLM based nuisance signal regression >
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

function [cleanData] = scrubData(data, mask, nuisance, probe, opts)
% Linear regression of input probes and nuisance regressors is performed on
% data in regions specified by the mask. Nuisance timeseries defined by
% nuisance regressors and co-efficients are summed and removed. For a
% similar function to determine residual signals after regression, see
% genGS.
%
% data: input timeseries data (i.e. 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% nuisance: and array of nuisance regressors (or a single regressor) having
% the same number of time-points as the input data.
%
% probe: an array of data-probes (explanatory variables) having the same
% number of time-points as the input data
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.figdir, opts.headers.ts, opts.resultsdir

warning('off');
global opts;
tf = class(data); 

if isfield(opts,'niiwrite'); else; opts.niiwrite = 0; end 
if isfield(opts,'save_cleaned'); else; opts.save_cleaned = 0; end

if isempty(probe)
    probe = zeros(length(nuisance));
else
    test1 = probe(1,:); test2 = probe(:,1);
    if length(test1) > length(test2); probe = probe'; end; clear test1 test2
    limits = [0 size(probe,1)];
    probemap = colormap(flip(brewermap(size(probe,2),'Spectral')));
end
if isempty(nuisance)
    nuisance = zeros(length(probe));
else
    test1 = nuisance(1,:); test2 = nuisance(:,1);
    if length(test1) > length(test2); nuisance = nuisance'; end; clear test1 test2
    limits = [0 size(nuisance,1)];
    nuisancemap = colormap(flip(brewermap(size(nuisance,2),'Spectral')));
end

for kk=1:size(nuisance,2)
% demean regressors
nuisance(:,kk) = rescale(demeanData(nuisance(:,kk)),-1,1);
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
if size(nuisance,2) > 1
    subplot(5,1,2); hold on; for ii=1:size(nuisance,2); plot(rescale(nuisance(:,ii),-1,1), 'Color', nuisancemap(ii,:)); end; title('rescaled nuisance regressor(s)'); xlim(limits);
else
    subplot(5,1,2); hold on; for ii=1:size(nuisance,2); plot(rescale(nuisance(:,ii),-1,1), 'Color', 'r'); end; title('rescaled nuisance regressor(s)'); xlim(limits);
end

subplot(5,1,3); plot(meanTimeseries(data,mask),'k'); title('Original Data'); xlim(limits)
subplot(5,1,4); plot(nanmean(nuis_TS,2),'k'); title('nuisance Mean'); xlim(limits);
%subplot(5,1,5); plot(meanTimeseries(cleanData,mask),'k'); title('Clean Data'); xlabel('datapoints'); xlim(limits);
subplot(5,1,5); title('Clean Data'); xlabel('datapoints'); xlim(limits);
plot(meanTimeseries(data,mask), 'k');
hold on
plot(meanTimeseries(cleanData,mask), 'r');
legend('before scrubbing', 'after scrubbing')
hold off


if isfield(opts,'figdir')
    saveas(gcf,fullfile(opts.figdir,'scrubData.fig'));
else
   saveas(gcf,fullfile(pwd,'scrubData.fig'));
end


if opts.save_cleaned
    if opts.niiwrite
        if isfield(opts.info,'rts'); else; opts.info.rts = opts.info.ts; end
        cd(opts.resultsdir);
        niftiwrite(cleanData,'cleanBOLD',opts.info.rts);
    else
        saveImageData(cleanData, opts.headers.ts, opts.resultsdir, 'cleanBOLD.nii.gz', 64);
    end
end
cleanData = cast(cleanData, tf);

end

