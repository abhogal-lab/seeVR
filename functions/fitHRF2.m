% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <fitHRF: Generates a series of HRFs and performs a best-fit GLM to input data >
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

if ispc
    if isfield(opts,'resultsdir'); else; opts.resultsdir = [pwd,'\']; end
else
    if isfield(opts,'resultsdir'); else; opts.resultsdir = [pwd,'/']; end
end
if ispc
    opts.hrfdir = [opts.resultsdir,'HRF\']; if exist(opts.hrfdir) == 7; else; mkdir(opts.hrfdir); end
else
    opts.hrfdir = [opts.resultsdir,'HRF/']; if exist(opts.hrfdir) == 7; else; mkdir(opts.hrfdir); end
end

t = cputime;

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

r_vec = zeros([1 length(coordinates)]);
HRF_vec = zeros([1 length(coordinates)]);
r2_map = zeros([1 x*y*z]);
HRF_map = zeros([1 x*y*z]);
onset_vec = zeros([1 length(coordinates)]);
disp_vec = zeros([1 length(coordinates)]);
under_vec = zeros([1 length(coordinates)]);
onset_map = zeros([1 x*y*z]);
disp_map = zeros([1 x*y*z]);
under_map = zeros([1 x*y*z]);

disp('Fitting HRF functions')
L = rescale(LegendreN(1,1:1:size(voxel_ts,2)));
parfor ii=1:size(voxel_ts,1)
    A = voxel_ts(ii,:);
    % with linear term
    if opts.include_linear
        C = [ones([length(A) 1]) L' A'];
        X = HRF_probe;
        regr_coef = C\HRF_probe';
        Y = regr_coef(3,:).*repmat(A,[ size(X,1) 1])' + regr_coef(2,:).*L' + regr_coef(1,:);        
    else
        C = [ones([length(A) 1])  A'];
        X = HRF_probe;
        regr_coef = C\HRF_probe';
        Y = regr_coef(2,:).*repmat(A,[ size(X,1) 1])' + regr_coef(1,:);        
    end
   

    SSE = zeros([1 size(X,1)]); SST = zeros([1 size(X,1)]);
    for jj=1:size(X,1)
        SSE(1,jj) = (norm(X(jj,:) - Y(:,jj)'))^2;
        SST(1,jj) = (norm(X(jj,:)- mean(X(jj,:))))^2;
    end
    
    R2 = 1 - SSE./SST;
    [M, I] = max(R2);
    r2_vec(1,ii) = M;
    HRF_vec(1,ii) = I;
    
    onset_vec(1,ii) = HRFidx(I,1);
    disp_vec(1,ii) = HRFidx(I,2);
    under_vec(1,ii) = HRFidx(I,3); 
end


r2_map(1, coordinates) = r2_vec; r2_map = reshape(r2_map, [x y z]);
HRF_map(1, coordinates) = HRF_vec; HRF_map = reshape(HRF_map, [x y z]);

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
if length(opts.onset) > 1; onset_map(1, coordinates) = onset_vec; onset_map = reshape(onset_map, [x y z]);
    if opts.niiwrite
        cd(opts.hrfdir);
        niftiwrite(onset_map,onsetName,opts.info.map);
    else
        saveImageData(onset_map,opts.headers.map,opts.hrfdir,onsetName,64);
    end
end
if length(opts.disp) > 1; disp_map(1, coordinates) = disp_vec; disp_map = reshape(disp_map, [x y z]);
    if opts.niiwrite
        cd(opts.hrfdir);
        niftiwrite(disp_map,dispName,opts.info.map);
    else
        saveImageData(disp_map,opts.headers.map,opts.hrfdir,dispName,64);
    end
end
if length(opts.under) > 1; under_map(1, coordinates) = under_vec; under_map = reshape(under_map, [x y z]);
    if opts.niiwrite
        cd(opts.hrfdir);
        niftiwrite(under_map,underName,opts.info.map);
    else
        saveImageData(under_map,opts.headers.map,opts.hrfdir,underName,64);
    end
end

if opts.niiwrite
    cd(opts.hrfdir);
    if opts.expHRF
        niftiwrite(HRF_map,'EXP_HRF_map',opts.info.map);
        niftiwrite(r2_map,'EXP_HRF_r2_map',opts.info.map);
    else
        niftiwrite(HRF_map,'GAMMA_HRF_map',opts.info.map);
        niftiwrite(r2_map,'GAMMA_HRF_r2_map',opts.info.map);
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

end