% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <fitTau: fits a convolved hemodynamic response function to extract dynamic (i.e. TAU) CVR information >
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
% this function is inspited by the manuscript from Poublanc et al.:
% Measuring cerebrovascular reactivity: the dynamic response to a step hypercapnic stimulus
% DOI: 10.1038/jcbfm.2015.114
%
% data: input timeseries data (i.e. 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% probe: usually PetCO2 or optimized regressor derived from lagCVR
% function. The proble is the 'reference' response
%
function [maps] = fitTau(probe, data, mask, opts)
global opts
tic

ty = class(data);
data = double(data);
mask = logical(mask);
if iscolumn(probe); probe = probe'; end

if isfield(opts,'interp_factor'); else; opts.interp_factor = 1; end        %factor by which to temporally interpolate data. Better for picking up lags between TR
if isfield(opts,'refine_tau'); else; opts.refine_tau = 1; end              %default refinement step for improved tau mapping using local information
if isfield(opts,'passes'); else; opts.passes = 10; end                     %maximum number of refinement passes
if isfield(opts,'win_size'); else; opts.win_size = 1; end                  %the number of voxels to consider around the voxel of interest when opts.refine_tau = 1;
if isfield(opts,'max_tau'); else; opts.max_tau = 300; end                  %maximum exponential dispersion time constant - data dependent
if isfield(opts,'save_unrefined'); else; opts.save_unrefined = 0; end          %save maps before refinement step to check effect

opts.dynamicdir = fullfile(opts.resultsdir,'tau'); mkdir(opts.dynamicdir);
[xx yy zz dyn] = size(data);

[voxel_ts, coordinates] = grabTimeseries(data, mask);

ts = zeros([length(coordinates),opts.interp_factor*dyn]);

if size(data,4) ~= length(probe) && opts.interp_factor > 1
    disp('interpolating data')
    t = opts.TR/opts.interp_factor:opts.TR/opts.interp_factor:size(data,4)*opts.TR;
    parfor ii = 1:length(coordinates)
        ts(ii,:) = interp(voxel_ts(ii,:), opts.interp_factor);
    end
    clear voxel_ts;
elseif size(data,4) == length(probe) && opts.interp_factor > 1
    disp('interpolating data and probe')
    t = opts.TR/opts.interp_factor:opts.TR/opts.interp_factor:size(data,4)*opts.TR;
    parfor ii = 1:length(coordinates)
        ts(ii,:) = interp(voxel_ts(ii,:), opts.interp_factor);
    end
    probe = interp(probe,opts.interp_factor);
    clear voxel_ts;
else
    t = opts.TR:opts.TR:size(data,4)*opts.TR;
    ts = voxel_ts; clear voxel_ts;
end

%rescale probe to help fit
input_probe = rescale(probe);

options = optimoptions('lsqcurvefit','Display','none','FunctionTolerance',1.0000e-8,...
    'StepTolerance', 1.0000e-8, 'MaxIter',150);

nr_params = 3;
b = nan([length(coordinates), nr_params]);

model = (@(a,t) a(1)*rescale(real(ifft(ifftshift(fftinput(input_probe).*fftexponential(a(2),t)))))+a(nr_params));

% bounds
lb = [ -20  0 -Inf];
ub = [ 30 opts.max_tau Inf];

parfor ii=1:length(coordinates)
    try
        b(ii,:) = lsqcurvefit(model,[range(ts(ii,:)) 20 0],t,ts(ii,:),lb,ub,options);
    catch
        printf(['error voxel ',int2str(ii)])
    end
end

b1_vec = zeros([1 xx*yy*zz]);
b2_vec = b1_vec; b3_vec = b1_vec;

b1_vec(1, coordinates) = b(:,1);
b2_vec(1, coordinates) = b(:,2);
b3_vec(1, coordinates) = b(:,3);

% refined analysis for clipped tau values

if opts.save_unrefined && opts.refine_tau
    
    b1_map = reshape(b1_vec, [xx yy zz]);
    b2_map = reshape(b2_vec, [xx yy zz]);
    b3_map = reshape(b3_vec, [xx yy zz]);
    
    if opts.niiwrite
        cd(opts.dynamicdir);
        niftiwrite(cast(mask.*b1_map,opts.mapDatatype),'exp_scaling_unrefined',opts.info.map);
        niftiwrite(cast(mask.*b2_map,opts.mapDatatype),'exp_tau_unrefined',opts.info.map);
        niftiwrite(cast(mask.*b3_map,opts.mapDatatype),'exp_offset_unrefined',opts.info.map);
    else
        saveImageData(mask.*b1_map, opts.headers.map, opts.dynamicdir, 'exp_scaling.nii.gz', 64);
        saveImageData(mask.*b2_map, opts.headers.map, opts.dynamicdir, 'exp_tau.nii.gz', 64);
        saveImageData(mask.*b3_map, opts.headers.map, opts.dynamicdir, 'exp_offset.nii.gz', 64);
    end
end

if opts.refine_tau
    
    b2_map = reshape(b2_vec, [xx yy zz]);
    iter = 1;
    passes = 0;
    max_tau = max(ceil(b2_map(:)));
    
    while iter
        passes = passes + 1
        max_map = zeros(size(mask));
        max_map(ceil(b2_map) == max_tau) = 1;
        
        %get new coordinates
        [~, newcoordinates] = grabTimeseries(data, max_map);
        clear b
        b = nan([length(newcoordinates), nr_params]);
        i = find(max_map); % find nonzero values in M
        [X,Y,Z] = ind2sub(size(max_map), i); clear i
        newTS = zeros([length(X) size(data,4)*opts.interp_factor]);
        
        for kk=1:length(X)
            try
                %extract timeseries and neighboring timeseries
                X_rng = X(kk)-opts.win_size:X(kk)+opts.win_size;
                Y_rng = Y(kk)-opts.win_size:Y(kk)+opts.win_size;
                Z_rng = Z(kk)-opts.win_size:Z(kk)+opts.win_size;
                %remove possiblity to have slices outside FOV
                X_rng(X_rng > size(data,1)) = []; X_rng(X_rng < 1) = [];
                Y_rng(Y_rng > size(data,2)) = []; Y_rng(Y_rng < 1) = [];
                Z_rng(Z_rng > size(data,3)) = []; Z_rng(Z_rng < 1) = [];
                
                %extract timeseries
                tmp =  data(X_rng,Y_rng,Z_rng,:);
                tmp = reshape(tmp, [size(tmp,1)*size(tmp,2)*size(tmp,3), size(tmp,4)]);
                tmp(find(tmp(:,1) == 0),:) = [];
                data(X(kk),Y(kk),Z(kk),:) = nanmean(tmp); %update data matrix for next iteration
                newTS(kk,:) = interp(nanmean(tmp),opts.interp_factor);
            catch
                disp('error; skipping voxel')
            end
        end
        
        clear tmp
        %fit for tau
        parfor ii=1:length(newcoordinates)
            try
                b(ii,:) = lsqcurvefit(model,[range(ts(ii,:)) 20 0],t,newTS(ii,:),lb,ub,options);
            catch
                disp(['error voxel ',int2str(ii)])
            end
        end
        b1_vec(1, newcoordinates) = b(:,1);
        b2_vec(1, newcoordinates) = b(:,2);
        b3_vec(1, newcoordinates) = b(:,3);
        
        %make new tau map for check
        b2_map = reshape(b2_vec, [xx yy zz]);
        max_map = zeros(size(b2_map));
        max_map(ceil(b2_map) == max_tau) = 1;
        perc = 100*(nnz(max_map)/nnz(mask));
        disp([int2str(perc), ' percent of voxels have clipped tau values'])
        if perc > 2
            if passes < opts.passes
                continue;
            else; iter = 0;
                disp('exceeded the maximum allowed passes')
                disp('... to increase passes set opts.passes to a higher value')
            end
        else
            iter = 0;
        end
    end
    
    b1_map = reshape(b1_vec, [xx yy zz]);
    b2_map = reshape(b2_vec, [xx yy zz]);
    b3_map = reshape(b3_vec, [xx yy zz]);
end
%% save maps and calculate stats

if opts.niiwrite
    cd(opts.dynamicdir);
    niftiwrite(cast(mask.*b1_map,opts.mapDatatype),'exp_scaling',opts.info.map);
    niftiwrite(cast(mask.*b2_map,opts.mapDatatype),'exp_tau',opts.info.map);
    niftiwrite(cast(mask.*b3_map,opts.mapDatatype),'exp_offset',opts.info.map);
else
    saveImageData(mask.*b1_map, opts.headers.map, opts.dynamicdir, 'exp_scaling.nii.gz', 64);
    saveImageData(mask.*b2_map, opts.headers.map, opts.dynamicdir, 'exp_tau.nii.gz', 64);
    saveImageData(mask.*b3_map, opts.headers.map, opts.dynamicdir, 'exp_offset.nii.gz', 64);
end

responseFits = zeros(length(coordinates), length(t));
b = [];
b(1,:) = b1_map(coordinates);
b(2,:) = b2_map(coordinates);
b(3,:) = b3_map(coordinates);
SSE = zeros(1,length(coordinates));
SST = zeros(1,length(coordinates));
r_corr = zeros(1,length(coordinates));

parfor ii=1:length(coordinates)
    responseFits(ii,:) = b(1,ii)*rescale(real(ifft(ifftshift(fftinput(input_probe).*fftexponential(b(2,ii),t)))))+b(3,ii);
    SSE(1,ii) = (norm(ts(ii,:) - responseFits(ii,:)))^2;
    SST(1,ii) = (norm(ts(ii,:)-mean(ts(ii,:))))^2;
    r_corr(1,ii) = corr(ts(ii,:)', responseFits(ii,:)');
end

R2 = 1 - SSE./SST; R2(R2 > 0) = 0;
cR2 = zeros(1, numel(mask)); r = zeros(1, numel(mask));
cR2(1, coordinates) = R2; r(1, coordinates) = r_corr;
cR2 = reshape(cR2, size(mask)); r = reshape(r, size(mask));

if opts.niiwrite
    niftiwrite(cast(mask.*cR2,opts.mapDatatype),'R2_map',opts.info.map);
    niftiwrite(cast(mask.*r,opts.mapDatatype),'r_map',opts.info.map);
    niftiwrite(cast(mask.*(r.^2),opts.mapDatatype),'expVariance_r2_map',opts.info.map);
else
    saveImageData(mask.*cR2, opts.headers.map, opts.corrCVRdir,'R2_map.nii.gz', 64);
    saveImageData(mask.*r, opts.headers.map, opts.corrCVRdir,'r_map.nii.gz', 64);
    saveImageData(mask.*(r.^2), opts.headers.map, opts.corrCVRdir,'expVariance_r2_map', 64);
end

maps.expHRF.scale = b1_map;
maps.expHRF.tau = b2_map;
maps.expHRF.offset = b3_map;
maps.expHRF.R2 = cR2;
maps.expHRF.r = r;
maps.expHRF.explainedVariance = r.^2;

toc

end