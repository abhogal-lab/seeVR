% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <denoiseData: performs wavelet-based denoising or temporal smoothing >
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
% Depending on the available toolboxes, this function will perform either a
% wavelet based de-noising or a temporal smoothing of timeseries data
%
% data: input timeseries data (i.e. 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% opts: options structure containing required variables for this specific
% function; i.e. for wavelets: opts.wdlevel, opts.family, opts.figdir. 
% for temporal smoothing: opts.method, opts.winsize (see MATLAB smoothdata
% function)
%
% denData: a denoised/smoothed version of the input data
function [denData] = denoiseData(data,mask,opts)   
global opts;
if isfield(opts,'wavelet'); else; opts.wavelet = 1; end
[xx yy zz dyn] = size(data);
tf = class(data);
if opts.wavelet
    try
if license('test', 'wavelet_toolbox')
    if isfield(opts,'wdlevel'); else; opts.wdlevel = 3; end
    if isfield(opts,'family'); else; opts.family = 'db4'; end
    [denData]  = wavDenoise(data, mask, opts);
end
    catch
    disp('wavelet toolbox unavalable; performing moving window mean smoothing')
    [voxel_ts, coordinates] = grabTimeseries(data, mask);
    if isfield(opts,'method'); else; opts.method = 'movmean'; end
    if isfield(opts,'winsize'); else; opts.winsize = round(0.05*size(data,4))'; end %5 percent window
    voxel_ts = smoothdata(voxel_ts,2,opts.method);
    denData = reshape(data, [xx*yy*zz dyn]);
    denData(coordinates,:) = voxel_ts;
    denData = reshape(denData, size(data));
    end
else
    [voxel_ts, coordinates] = grabTimeseries(data, mask);
    if isfield(opts,'method'); else; opts.method = 'movmean'; end
    if isfield(opts,'winsize'); else; opts.winsize = round(0.05*size(data,4))'; end %5 percent window
    voxel_ts = smoothdata(voxel_ts,2,opts.method);
    denData = reshape(data, [xx*yy*zz dyn]);
    denData(coordinates,:) = voxel_ts;
    denData = reshape(denData, size(data));  
end
if opts.verbose

figure; hold on;
xdata = [opts.TR:opts.TR:opts.TR*size(data,4)];
plot(xdata, meanTimeseries(data,mask), 'k', 'LineWidth',2);
plot(xdata, meanTimeseries(denData,mask), 'c', 'LineWidth',2);
title('noise removal');
ylabel('absolute signal'); xlabel('time (s)'); legend('original', 'de-noised')
set(gcf, 'Units', 'pixels', 'Position', [200, 500, 600, 160]); xlim([0 xdata(end)])
hold off
    if isfield(opts, 'figdir')
        saveas(gcf,[opts.figdir,'wavDenoise_data.fig']);
    else
        if ispc
            saveas(gcf,[pwd,'\','wavDenoise_data.fig']);
        else
            saveas(gcf,[pwd,'/','wavDenoise_data.fig']);
        end
    end
end
denData = cast(denData,tf);
end

