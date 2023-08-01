% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <fitDynamics: fits a convolved hemodynamic response function to extract dynamic (i.e. TAU) CVR information >
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
% probe: usually PetCO2 or optimized regressor derived from lagCVR
% function. The proble is the 'reference' response
%
% nuisance: set of nuisance parameters for regression

function [maps] = fitTau(probe, data, mask, opts)
global opts
ty = class(data);
data = double(data);
mask = logical(mask);
if iscolumn(probe); probe = probe'; end

if isfield(opts,'interp_factor'); else; opts.interp_factor = 2; end       %factor by which to temporally interpolate data. Better for picking up lags between TR
if isfield(opts,'maxTau'); else; opts.maxTau = 300; end                   %maximum exponential dispersion time constant - data dependent

opts.dynamicdir = fullfile(opts.resultsdir,'tau'); mkdir(opts.dynamicdir);
[xx yy zz dyn] = size(data);

[voxel_ts, coordinates] = grabTimeseries(data, mask);

ts = zeros([length(coordinates),opts.interp_factor*dyn]);

if size(data,4) ~= length(probe) && opts.interp_factor > 1
    disp('interpolating data')
    t = opts.TR/opts.interp_factor:opts.TR/opts.interp_factor:size(data,4)*opts.TR;
    parfor ii = 1:length(coordinates)
        ts(ii,:) = interp(voxel_ts(ii,:),opts.interp_factor);
    end
    clear voxel_ts;
elseif size(data,4) == length(probe) && opts.interp_factor > 1
    disp('interpolating data and probe')
    t = opts.TR/opts.interp_factor:opts.TR/opts.interp_factor:size(data,4)*opts.TR;
    parfor ii = 1:length(coordinates)
        ts(ii,:) = interp(voxel_ts(ii,:),opts.interp_factor);
    end
    probe = interp(probe,opts.interp_factor);
    clear voxel_ts;
else
    t = opts.TR:opts.TR:size(data,4)*opts.TR;
    ts = voxel_ts; clear voxel_ts;
end

options = optimoptions('lsqcurvefit','Display','none','FunctionTolerance',1.0000e-8,...
    'StepTolerance', 1.0000e-8, 'MaxIter',150);

nr_params = 3

b = nan([length(coordinates), nr_params]);

model = (@(a,t) a(1)*rescale(real(ifft(ifftshift(fftinput(probe).*fftexponential(a(2),t)))))+a(nr_params));
a0 = [0 20 0];
lb = [ 0  0 -10];
ub = [ 30 opts.maxTau 20];

tic

parfor ii=1:length(coordinates)
    try
        b(ii,:) = lsqcurvefit(model,[range(ts(ii,:)) 20 0],t,ts(ii,:),lb,ub,options);
    catch
        %continue
        printf(['error voxel ',int2str(ii)])
    end
    %         y_fit = model(b(ii,:),t);
    %         % Plot the results
    %         figure
    %         plot(t,(ts(ii,:)),'r-',t,(y_fit),'b-');
    %         legend('Observed Data','Fitted Model');
    %         pause(0.1)
end

b1_vec = zeros([1 xx*yy*zz]);
b2_vec = b1_vec; b3_vec = b1_vec;

b1_vec(1, coordinates) = b(:,1); b1_map = reshape(b1_vec, [xx yy zz]);
b2_vec(1, coordinates) = b(:,2); b2_map = reshape(b2_vec, [xx yy zz]);
b3_vec(1, coordinates) = b(:,3); b3_map = reshape(b3_vec, [xx yy zz]);

if opts.niiwrite
    cd(opts.dynamicdir);
    niftiwrite(cast(mask.*b1_map,opts.mapDatatype),'exp_scaling',opts.info.map);
    niftiwrite(cast(mask.*b2_map,opts.mapDatatype),'exp_tau',opts.info.map);
    niftiwrite(cast(mask.*b3_map,opts.mapDatatype),'exp_offset',opts.info.map);
else
    saveImageData(mask.*b1_map, opts.headers.map, opts.dynamicdir, 'exp_scaling.nii.gz', datatype);
    saveImageData(mask.*b2_map, opts.headers.map, opts.dynamicdir, 'exp_tau.nii.gz', datatype);
    saveImageData(mask.*b3_map, opts.headers.map, opts.dynamicdir, 'exp_offset.nii.gz', datatype);
    
end

maps.expHRF.scale = b1_map;
maps.expHRF.tau = b2_map;
maps.expHRF.offset = b3_map;

end