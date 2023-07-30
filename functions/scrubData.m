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
% the same number of time-points as the input data. Generally these will be
% the rotations and translations derived from motion correction
%
% probe: an array of data-probes (explanatory variables) having the same
% number of time-points as the input data
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.figdir, opts.headers.ts, opts.resultsdir

warning('off');
global opts;
tf = class(data);
data = double(data);
nuisance = double(nuisance);
probe = double(rescale(probe));

[x,y,z,dyn] = size(data);
[voxel_ts, coordinates] = grabTimeseries(data, mask);

%function options
if isfield(opts,'niiwrite'); else; opts.niiwrite = 0; end
if isfield(opts,'save_cleaned'); else; opts.save_cleaned = 0; end
if isfield(opts,'disperse_probe'); else; opts.disperse_probe = 1; end %convolves the input probe with a few exponentials to create additional EVs

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
end


%prepare nuisance probes
[np0] = prepNuisance(nuisance,probe, opts);
nuisancemap = colormap(flip(brewermap(size(np0,2),'Spectral')));

if opts.disperse_probe
opts.disp = [5 15 25 35]; %dispersion
[~,~,probe] = convHRF(probe, opts);
probe = probe';
end

probemap = colormap(flip(brewermap(size(probe,2),'Spectral')));
probesize = size(probe,2);

%perform regression
np = [probe np0];

D = [ones([length(np) 1]) np];
np_coef = D\voxel_ts';

%nnintcp = squeeze(np_coef(1,:));
nuis = squeeze(np_coef(2+probesize:end,:));
nuis_TS = zeros(size(voxel_ts))';

parfor ii=1:size(np0,2)
    tmp = np0(:,ii)*nuis(ii,:);
    nuis_TS = nuis_TS+tmp;
end

switch ndims(data)
    case 4
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

%plot results
figure;
if size(probe,2) > 1
    subplot(5,1,1); hold on; for ii=1:size(probe,2); plot(probe(:,ii), 'Color', probemap(ii,:)); end; title('data probe(s)'); xlim(limits);
else
    subplot(5,1,1); hold on; for ii=1:size(probe,2); plot(probe(:,ii), 'Color', 'k'); end; title('data probe(s)'); xlim(limits);
end
if size(np0,2) > 1
    subplot(5,1,2); hold on; for ii=1:size(np0,2); plot(np0(:,ii), 'Color', nuisancemap(ii,:)); end; title('nuisance regressor(s)'); xlim(limits);
else
    subplot(5,1,2); hold on; for ii=1:size(np0,2); plot(np0(:,ii), 'Color', 'r'); end; title('nuisance regressor(s)'); xlim(limits);
end

subplot(5,1,3); plot(meanTimeseries(data,mask),'k'); title('Original Data'); xlim(limits)
subplot(5,1,4); plot(nanmean(nuis_TS,2),'k'); title('nuisance Mean'); xlim(limits);
subplot(5,1,5); title('Clean Data'); xlabel('datapoints');
plot(meanTimeseries(data,mask), 'k');
hold on
plot(meanTimeseries(cleanData,mask), 'r'); xlim(limits);
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

