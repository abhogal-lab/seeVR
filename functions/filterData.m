% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <filterData: performs a variety of smoothing operations depending on parameter inputs >
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
% 1) 2D/3D gaussian smoothing
% 2) 2D/3D guideimd bilateral smoothing for edge preservation (default option)
% 3) TIPS (time–intensity profile similarity)smoothing
% (https://iopscience.iop.org/article/10.1088/0031-9155/56/13/008)
%
%   guideim: should be a gray-scale/single channel image with structural
%   information (i.e. the first volume in the original BOLD timeseries)
%
%   data: should be a gray-scale/single channel image of data to be
%   smoothed (i.e. baseline-normalized BOLD data - %BOLD)
%
%   opts.filter: the selected filtering method - 'tips', 'bilateral',
%   'gaussian'
%
%   opts.FWHM = FWHM of smoothing kernel in mm
%
%   opts.sigma_range: parameter that defines smoothing threshold - data
%   dependent
%
% *************************************************************************

function [filtered_data] = filterData(data,guideim,mask,opts)
data(isnan(data)) = 0;
warning('off');
global opts;
type = class(data);
if isfield(opts,'filter'); else; opts.filter = 'bilateral'; end
if isfield(opts,'spatialdim'); else; opts.spatialdim = 2; end
if isfield(opts,'FWHM'); else; opts.FWHM = [4 4 4]; end
if isfield(opts,'sigma_range') 
    
%do nothing
else 
    if ndims(data) > 3
    opts.sigma_range = 10*ROIstd(mean(data(:,:,:,1:10),4),mask); 
    else
        tmp = data.*mask;
        tmp = tmp(:);
        tmp(tmp == 0) = [];
    opts.sigma_range = 10*std(tmp); 
    end
end
%FWHM = voxelSize*sigma_spatial*2.355;
%sigma_spatial = FWHM/(2.355*voxelsize)

if numel(opts.FWHM) == 1
    tmp = opts.FWHM;
    switch opts.spatialdim
        case 2        
        opts.FWHM = [tmp tmp]; 
        case 3
        opts.FWHM = [tmp tmp tmp];
    end
end

if opts.spatialdim == 2
    opts.sigma_spatial = [opts.FWHM(1)/(opts.voxelsize(1)*2.355) opts.FWHM(2)/(opts.voxelsize(2)*2.355)  0];
else
    opts.sigma_spatial = [opts.FWHM(1)/(opts.voxelsize(1)*2.355) opts.FWHM(2)/(opts.voxelsize(2)*2.355) opts.FWHM(3)/(opts.voxelsize(3)*2.355)];
end

filtered_data = zeros(size(data));

if isfield(opts,'sigma_range'); else
    %opts.sigma_range = sqrt(2)*ROIstd(data(:,:,:,1),mask); %just a guess
    opts.sigma_range = ROIstd(guideim,mask); %just a guess
end

%run specified filter

opts.sigma_spatial = double(opts.sigma_spatial);
opts.sigma_range = double(opts.sigma_range);

switch opts.filter
    
    case 'imguided'
        if isfield(opts,'nHood'); else; opts.nHood = [3 3]; end
        opts.spatialdim = 2;
        disp('This is a 2D implementation; opts.spatialdim updated to 2')
        
        for ii=1:size(data,4)
            filtered_data(:,:,:,ii) = imguidedfilter(squeeze(data(:,:,:,ii)), guideim, 'NeighborhoodSize', opts.nHood);
        end
        
    case 'bilateral'
        
        mask = uint8(mask);
        %setup CPUs
        n_cpu = java.lang.Runtime.getRuntime().availableProcessors();
        n_threads = max(n_cpu - 1, 1);
        setenv('OMP_NUM_THREADS', num2str(n_threads));
        
        tdata = single(permute(data,[4 1 2 3]));
        
        if ispc
            filtered_data = winbilatFilter(tdata, single(guideim), double(opts.sigma_spatial), opts.sigma_range, mask);
        else
            if ismac
                filtered_data = macbilatFilter(tdata,single(guideim),opts.sigma_spatial, opts.sigma_range, mask);
            else
                filtered_data = bilatFilter(tdata,single(guideim),opts.sigma_spatial, opts.sigma_range, mask);
            end
        end
        filtered_data = permute(filtered_data, [2 3 4 1]);
        
    case 'tips'
        
        mask = uint8(mask);
        %setup CPUs
        n_cpu = java.lang.Runtime.getRuntime().availableProcessors();
        n_threads = max(n_cpu - 1, 1);
        setenv('OMP_NUM_THREADS', num2str(n_threads));
        
        tdata = single(permute(data,[4 1 2 3]));
        
        if ispc
            filtered_data = wintipsFilter(tdata,opts.sigma_spatial, opts.sigma_range, uint8(mask));
        else
            if ismac
                filtered_data = mactipsFilter(tdata,opts.sigma_spatial, opts.sigma_range, mask);
            else
                filtered_data = tipsFilter(tdata,opts.sigma_spatial, opts.sigma_range, mask);
            end
        end
        filtered_data = permute(filtered_data, [2 3 4 1]);
        
    case 'gaussian'
        
        %to emulate gaussian smoothing, set very high range
        opts.sigma_range = 1000000000;
        
        mask = uint8(mask);
        %setup CPUs
        n_cpu = java.lang.Runtime.getRuntime().availableProcessors();
        n_threads = max(n_cpu - 1, 1);
        setenv('OMP_NUM_THREADS', num2str(n_threads));
        
        tdata = single(permute(data,[4 1 2 3]));
        
        if ispc
            filtered_data = winbilatFilter(tdata,single(guideim),opts.sigma_spatial, opts.sigma_range, mask);
        else
            if ismac
                filtered_data = macbilatFilter(tdata,single(guideim),opts.sigma_spatial, opts.sigma_range, mask);
            else
                filtered_data = bilatFilter(tdata,single(guideim),opts.sigma_spatial, opts.sigma_range, mask);
            end
        end
        filtered_data = permute(filtered_data, [2 3 4 1]);
        
end
filtered_data = eval([type,'(filtered_data)']);
end
