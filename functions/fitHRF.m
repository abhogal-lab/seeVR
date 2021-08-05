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
function [maps] = fitHRF(mask,data,probe,opts)
global opts
warning('off')

if isfield(opts,'verbose'); else; opts.verbose = 0; end                     %turn on/off select command output
if isfield(opts,'prewhite'); else; opts.prewhite = 0; end                   %pre-whiten data
if isfield(opts,'interp_factor'); else; opts.interp_factor = 1; end         %interpolate timeseries
if isfield(opts,'expHRF'); else; opts.expHRF = 0; end                       %convolve with exponential HRF (default = 0; use double gamma)


if ispc
    if isfield(opts,'resultsdir'); else; opts.resultsdir = [pwd,'\']; end
else
    if isfield(opts,'resultsdir'); else; opts.resultsdir = [pwd,'/']; end
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
[~,HRF_probe] = convEXP(probe, opts); %Generate a exponential HRF 
else
if opts.verbose; disp('generating HRF maps based on double gamma'); opts.plot = 1; end
[~,HRF_probe] = convHRF(probe, opts); %Generate a double gamma HRF
end

r_vec = zeros([1 length(coordinates)]);
HRF_vec = zeros([1 length(coordinates)]);
r2_map = zeros([1 x*y*z]);
HRF_map = zeros([1 x*y*z]);

disp('Fitting HRF functions')

parfor ii=1:size(voxel_ts,1)
    A = voxel_ts(ii,:);
    % with linear term
    C = [ones([length(A) 1])  A'];
    regr_coef = C\HRF_probe';
    
    %original observations
    X = HRF_probe;
    %fitted values
    Y = regr_coef(2,:).*repmat(A,[ size(X,1) 1])' + regr_coef(1,:);
    
    SSE = zeros([1 size(X,1)]); SST = zeros([1 size(X,1)]);
    for jj=1:size(X,1)
        SSE(1,jj) = (norm(X(jj,:) - Y(:,jj)'))^2;
        SST(1,jj) = (norm(X(jj,:)- mean(X(jj,:))))^2;
    end
    
    R2 = 1 - SSE./SST;
    [M, I] = max(R2);
    r2_vec(1,ii) = M;
    HRF_vec(1,ii) = I;
end


r2_map(1, coordinates) = r2_vec; r2_map = reshape(r2_map, [ x y z]);
HRF_map(1, coordinates) = HRF_vec; HRF_map = reshape(HRF_map, [ x y z]);

if opts.expHRF
saveImageData(HRF_map,opts.headers.map,opts.resultsdir,'expHRF_map.nii.gz',64);
saveImageData(r2_map,opts.headers.map,opts.resultsdir,'expHRF_r2_map.nii.gz',64);    
else
if length(opts.disp) > length(opts.onset)
saveImageData(HRF_map,opts.headers.map,opts.resultsdir,'dispHRF_map.nii.gz',64);
saveImageData(r2_map,opts.headers.map,opts.resultsdir,'dispHRF_r2_map.nii.gz',64);
elseif length(opts.disp) < length(opts.onset)
saveImageData(HRF_map,opts.headers.map,opts.resultsdir,'onsetHRF_map.nii.gz',64);
saveImageData(r2_map,opts.headers.map,opts.resultsdir,'onsetHRF_r2_map.nii.gz',64);
else
saveImageData(HRF_map,opts.headers.map,opts.resultsdir,'HRF_map.nii.gz',64);
saveImageData(r2_map,opts.headers.map,opts.resultsdir,'HRF_r2_map.nii.gz',64);   
end
end

maps.HRF = HRF_map;
maps.r2 = r2_map;
if opts.verbose; disp(['HRF maps were created in: ',int2str(cputime-t),' seconds']); end

end