% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <genGS: GLM based approach to regress out nuisance and data signals in provide residual timeseries data >
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
% This function uses input nuissance regressors and data regressors
% (i.e. response models) to remove input data of all explainable
% signal responses. The result is a residual 'global' signal that can be
% added as a further nuisance regressor (see scrubData) or used as 'pseudo
% resting-state' data.
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
% function; i.e. opts.figdir
%
% cleanData: residual data where nuisance and probe regressor explained
% signal variance has been removed. This data retains the linear offset
% from the GLM regression; in this way genGS can also be used to clean data
% of unwanted signal contributions (see scrubData)
%
% res_ts: the mean timeseries (cleanData) signal calculated within the ROI defined by
% the input mask
%
% resData: residual data after nuisance and probe regressor explained
% signal variance has been removed. Here the linear offset is also removed.
function [cleanData,res_ts,resData] = genGS(data, mask, nuisance, probe, opts)

warning('off');
global opts

if isempty(probe); clear probe; end
if isempty(nuisance); clear nuisance; end

%check arguments
if exist('probe','var') == 1 && exist('nuisance','var') == 1
    test1 = probe(1,:); test2 = probe(:,1);
    if length(test1) > length(test2); probe = probe'; end; clear test1 test2
    test1 = nuisance(1,:); test2 = nuisance(:,1);
    if length(test1) > length(test2); nuisance = nuisance'; end; clear test1 test2
    disp('nuisance and data probes found');
    check = 1;
elseif exist('probe','var') == 1
    test1 = probe(1,:); test2 = probe(:,1);
    if length(test1) > length(test2); probe = probe'; end; clear test1 test2
    disp('nuisance and data probes found');
    check = 2;
elseif exist('nuisance','var') == 1
    test1 = nuisance(1,:); test2 = nuisance(:,1);
    if length(test1) > length(test2); nuisance = nuisance'; end; clear test1 test2
    disp('nuisance and data probes found');
    check = 3;
end

switch check
    case 1
        np = [nuisance probe];
    case 2
        np = probe;
    case 3
        np = nuisance;
end

[voxel_ts, coordinates] = grabTimeseries(data, mask);

D = [ones([length(np) 1]) np];
np_coef = []; nintcp = []; nuis = []; nuis_TS = []; ncombi_TS=[];
np_coef = D\voxel_ts';
nintcp = squeeze(np_coef(1,:));
nuis = squeeze(np_coef(2:end,:));

nuis_TS = zeros(size(voxel_ts))';
for ii=1:size(nuis,1)
    tmp = np(:,ii)*nuis(ii,:);
    nuis_TS = nuis_TS + tmp;
end

offset = ones([length(np) 1])*nintcp;
ncombi_TS = nuis_TS;

switch ndims(data)
    case 4
        [x,y,z,dyn] = size(data);
        res_data = voxel_ts - ncombi_TS';
        resData = zeros([x*y*z,dyn]);
        resData(coordinates,:) = res_data;
        resData = reshape(resData,[x y z dyn]);
        
        offBOLD = zeros([x*y*z,dyn]);
        offBOLD(coordinates,:) = offset';
        offBOLD = reshape(offBOLD,[x y z dyn]);
        cleanData = data - resData + offBOLD;
    case 3
        [x,y,dyn] = size(data);
        res_data = voxel_ts - ncombi_TS';
        resData = zeros([x*y,dyn]);
        resData(coordinates,:) = res_data;
        resData = reshape(resData,[x y dyn]);
        
        offBOLD = zeros([x*y,dyn]);
        offBOLD(coordinates,:) = offset';
        offBOLD = reshape(offBOLD,[x y dyn]);
        cleanData = data - resData + offBOLD;
end
%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%
switch check
    case 1
        figure;
        probemap = plasma(size(probe,2));
        nuisancemap = viridis(size(nuisance,2)); nuisancemap = flip(nuisancemap,1);
        limits = [0 size(probe,1)];
        subplot(5,1,1); plot(meanTimeseries(data,mask), 'k'); title('original data'); xlim(limits);
        subplot(5,1,2); hold on; for ii=1:size(probe,2); plot(probe(:,ii), 'Color', probemap(ii,:)); end; title('data probes'); xlim(limits);
        subplot(5,1,3); hold on; for ii=1:size(nuisance,2); plot(rescale(nuisance(:,ii),-1,1), 'Color', nuisancemap(ii,:)); end; title('rescaled nuisance regressors'); xlim(limits);
        subplot(5,1,4); plot(meanTimeseries(cleanData,mask), 'k'); title('explained data'); xlim(limits);
        subplot(5,1,5); plot(meanTimeseries(resData-offBOLD,mask), 'k'); title('residual data'); xlim(limits);
    case 2
        figure;
        probemap = plasma(size(probe,2));
        limits = [0 size(probe,1)];
        subplot(4,1,1); plot(meanTimeseries(data,mask), 'k'); title('original data'); xlim(limits);
        subplot(4,1,2); hold on; for ii=1:size(probe,2); plot(probe(:,ii), 'Color', probemap(ii,:)); end; title('data probes'); xlim(limits);
        subplot(4,1,3); plot(meanTimeseries(cleanData,mask), 'k'); title('explained data'); xlim(limits);
        subplot(4,1,4); plot(meanTimeseries(resData-offBOLD,mask), 'k'); title('residual data'); xlim(limits);
    case 3
        figure;
        nuisancemap = viridis(size(nuisance,2)); nuisancemap = flip(nuisancemap,1);
        limits = [0 size(nuisance,1)];
        subplot(4,1,1); plot(meanTimeseries(data,mask), 'k'); title('original data'); xlim(limits);
        subplot(4,1,2); hold on; for ii=1:size(nuisance,2); plot(rescale(nuisance(:,ii),-1,1), 'Color', nuisancemap(ii,:)); end; title('rescaled nuisance regressors'); xlim(limits);
        subplot(4,1,3); plot(meanTimeseries(cleanData,mask), 'k'); title('explained data'); xlim(limits);
        subplot(4,1,4); plot(meanTimeseries(resData-offBOLD,mask), 'k'); title('residual data'); xlim(limits);
end

resData = resData-offBOLD;

if isfield(opts,'figdir')
    saveas(gcf,fullfile(opts.figdir,'genGS.fig'));
else
    saveas(gcf,fullfile(pwd,'genGS.fig'));
end

res_ts = meanTimeseries(resData,mask);

end

