% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <fitHRF2: Generates a series of HRFs and performs a best-fit GLM to input data >
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
%
%
function [maps] = fitHRF2(mask,data,probe,opts)
global opts
warning('off')

if isfield(opts,'verbose'); else; opts.verbose = 0; end                    %turn on/off select command output
if isfield(opts,'prewhite'); else; opts.prewhite = 0; end                  %pre-whiten data
if isfield(opts,'interp_factor'); else; opts.interp_factor = 1; end        %interpolate timeseries
if isfield(opts,'expHRF'); else; opts.expHRF = 0; end                      %convolve with exponential HRF (default = 0; use double gamma)
if isfield(opts,'niiwrite'); else; opts.niiwrite = 0; end                  %depending on how data is loaded this can be set to 1 to use native load/save functions
if isfield(opts,'include_linear'); else; opts.include_linear = 0; end      %add a liner term to the HRF fit - this can possible account for strange behavior due to the convolution
   
if isfield(opts,'resultsdir'); else; opts.resultsdir = pwd; end

opts.hrfdir = fullfile(opts.resultsdir,'HRF'); if exist(opts.hrfdir, 'dir') == 7; else; mkdir(opts.hrfdir); end
[x, y, z, dyn] = size(data);


%% grab coordinates
[clean_voxel_ts, coordinates] = grabTimeseries(data, mask);

if opts.prewhite
    clean_voxel_ts = clean_voxel_ts'; parfor ii=1:size(coordinates,1); [clean_voxel_ts(:,ii), ~, ~] = prewhiten(clean_voxel_ts(:,ii)); end; clean_voxel_ts = clean_voxel_ts';
    [probe, ~, ~] = prewhiten(probe);
end

%interpolate data - not sure if needed for HRF fitting
if opts.interp_factor > 1
    voxel_ts=zeros([length(coordinates) opts.interp_factor*dyn]);
    if length(probe)/size(data,ndims(data)) == opts.interp_factor
        %do nothing since the input probe has likel already been interpolated
    elseif length(probe)/size(data,ndims(data)) == 1
        probe = interp(probe,opts.interp_factor);
    else
        error('check that the length of the input probe and data are the same; exiting function')
        return
    end
    parfor ii = 1:length(coordinates)
        voxel_ts(ii,:) = interp(clean_voxel_ts(ii,:),opts.interp_factor);
    end
else
    voxel_ts = clean_voxel_ts;
end
%% Generate HRF functions for fitting

hrf_map = zeros([1 x*y*z]); hrf_estimate_map = zeros([1 x*y*z]);

if opts.expHRF
    if opts.verbose; disp('generating HRF maps based on exponential'); opts.plot = 1; end
    opts.under = 1;
    [HRFidx,HRF_probe] = convEXP2(probe, opts); %Generate a exponential HRF
else
    if opts.verbose; disp('generating HRF maps based on double gamma'); opts.plot = 1; end
    [HRFidx, HRF_probe] = convHRF2(probe, opts); %Generate a double gamma HRF
end

r2_map = zeros([1 x*y*z]);
HRF_map = zeros([1 x*y*z]);
onset_map = zeros([1 x*y*z]);
disp_map = zeros([1 x*y*z]);
under_map = zeros([1 x*y*z]);
beta_vec = zeros([1 length(coordinates)]);
beta_map = zeros([1 x*y*z]);

disp('Fitting HRF functions')

if opts.include_linear
    L = rescale(legendreP(1,1:1:size(voxel_ts,2)),-1,1);
    regr_coef = zeros([size(HRF_probe,1) 3 length(coordinates)]);
    parfor ii=1:size(HRF_probe,1)
        A = HRF_probe(ii,:);
        C = [ones([length(A(1,~isnan(A))) 1]) L(1,~isnan(A))' A(1,~isnan(A))'];
        regr_coef(ii,:,:)= C\voxel_ts(:,~isnan(A))';
    end
    
    %reconstruct signals to calculate R2
    R2 = zeros([1 size(HRF_probe,1)]);
    maxindex = zeros([1 size(voxel_ts,1)]);
    beta = zeros([1 size(voxel_ts,1)]);
    rsquared = zeros([1 size(voxel_ts,1)]);
    for ii = 1:size(voxel_ts,1)
        tcoef = regr_coef(:,:,ii);
        X =  voxel_ts(ii,:);
        for jj=1:size(HRF_probe,1)
            tY = tcoef(jj,3).*HRF_probe(jj,:) + tcoef(jj,2).*L + tcoef(jj,1);
            tmp = corrcoef(tY,X);
            R2(jj) = tmp(1,2)^2;
        end
        [M, I] = max(R2');
        maxindex(1,ii) = I;
        rsquared(1,ii) = M;
        beta(1,ii) = regr_coef(I,2,ii);
    end
else %do not add a linear term to the regression
    
    regr_coef = zeros([size(HRF_probe,1) 2 length(coordinates)]);
    parfor ii=1:size(HRF_probe,1)
        A = HRF_probe(ii,:);
        C = [ones([length(A(1,~isnan(A))) 1]) A(1,~isnan(A))'];
        regr_coef(ii,:,:)= C\voxel_ts(:,~isnan(A))';
    end
    
    %reconstruct signals to calculate R2
    R2 = zeros([1 size(HRF_probe,1)]);
    maxindex = zeros([1 size(voxel_ts,1)]);
    beta = zeros([1 size(voxel_ts,1)]);
    rsquared = zeros([1 size(voxel_ts,1)]);
    for ii = 1:size(voxel_ts,1)
        tcoef = regr_coef(:,:,ii);
        X =  voxel_ts(ii,:);
        for jj=1:size(HRF_probe,1)
            tY = tcoef(jj,2).*HRF_probe(jj,:) + tcoef(jj,1);
            tmp = corrcoef(tY,X);
            R2(jj) = tmp(1,2)^2;
        end
        [M, I] = max(R2');
        maxindex(1,ii) = I;
        rsquared(1,ii) = M;
        beta(1,ii) = regr_coef(I,2,ii);
    end
end

onset_vec(1,:) = HRFidx(maxindex,1);
disp_vec(1,:) = HRFidx(maxindex,2);
under_vec(1,:) = HRFidx(maxindex,3);

beta_map(1, coordinates) = beta; beta_map = reshape(beta_map, [x y z]);
r2_map(1, coordinates) = rsquared; r2_map = reshape(r2_map, [x y z]);
HRF_map(1, coordinates) = maxindex; HRF_map = reshape(HRF_map, [x y z]);

if opts.expHRF
    if opts.niiwrite
        onsetName = 'EXP_onset_map';
        dispName = 'EXP_dispersion_map';
        underName = 'EXP_undershoot_map';
    else
        onsetName = 'EXP_onset_map.nii';
        dispName = 'EXP_dispersion_map.nii';
        underName = 'EXP_undershoot_map.nii';
    end
else
    if opts.niiwrite
        onsetName = 'GAMMA_onset_map';
        dispName = 'GAMMA_dispersion_map';
        underName = 'GAMMA_undershoot_map';
    else
        onsetName = 'GAMMA_onset_map.nii';
        dispName = 'GAMMA_dispersion_map.nii';
        underName = 'GAMMA_undershoot_map.nii';
    end
end
if length(opts.onset) > 1; onset_map(coordinates) = onset_vec; onset_map = reshape(onset_map, [x y z]);
    if opts.niiwrite
        cd(opts.hrfdir);
        niftiwrite(cast(onset_map,opts.mapDatatype),onsetName,opts.info.map);
    else
        saveImageData(onset_map,opts.headers.map,opts.hrfdir,onsetName,64);
    end
end
if length(opts.disp) > 1; disp_map(coordinates) = disp_vec; disp_map = reshape(disp_map, [x y z]);
    if opts.niiwrite
        cd(opts.hrfdir);
        niftiwrite(cast(disp_map,opts.mapDatatype),dispName,opts.info.map);
    else
        saveImageData(disp_map,opts.headers.map,opts.hrfdir,dispName,64);
    end
end
if length(opts.under) > 1; under_map(coordinates) = under_vec; under_map = reshape(under_map, [x y z]);
    if opts.niiwrite
        cd(opts.hrfdir);
        niftiwrite(cast(under_map,opts.mapDatatype),underName,opts.info.map);
    else
        saveImageData(under_map,opts.headers.map,opts.hrfdir,underName,64);
    end
end

if opts.niiwrite
    cd(opts.hrfdir);
    if opts.expHRF
        niftiwrite(cast(HRF_map,opts.mapDatatype),'EXP_HRF_map',opts.info.map);
        niftiwrite(cast(r2_map,opts.mapDatatype),'EXP_HRF_r2_map',opts.info.map);
        niftiwrite(cast(beta_map,opts.mapDatatype),'EXP_HRF_beta_map',opts.info.map);        
    else
        niftiwrite(cast(HRF_map,opts.mapDatatype),'GAMMA_HRF_map',opts.info.map);
        niftiwrite(cast(r2_map,opts.mapDatatype),'GAMMA_HRF_r2_map',opts.info.map);
        niftiwrite(cast(beta_map,opts.mapDatatype),'GAMMA_HRF_beta_map',opts.info.map);       
    end
else
    if opts.expHRF
        saveImageData(HRF_map,opts.headers.map,opts.hrfdir,'EXP_HRF_map.nii.gz',64);
        saveImageData(r2_map,opts.headers.map,opts.hrfdir,'EXP_HRF_r2_map.nii.gz',64);
    else
        saveImageData(HRF_map,opts.headers.map,opts.hrfdir,'GAMMA_HRF_map.nii.gz',64);
        saveImageData(r2_map,opts.headers.map,opts.hrfdir,'GAMMA_HRF_r2_map.nii.gz',64);
    end
end

maps.HRF = HRF_map;
maps.r2 = r2_map;
maps.onset = onset_map;
maps.disp = disp_map;
maps.under = under_map;
maps.beta = beta_map;

save(fullfile(opts.resultsdir,'HRF_maps.mat'), 'maps');
end