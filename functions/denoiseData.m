function [denData] = denoiseData(data,mask,opts)
%written by Alex Bhogal a.bhogal@umcutrecht.nl
%This functions performs a wavelet denoising on timeseries data. If wavelet
%toolbox is unavailable then a standard smoothing is performed based on
%used specified parameters

if license('test', 'wavelet_toolbox')
    if isfield(opts,'wdlevel'); else; opts.wdlevel = 3; end
    if isfield(opts,'family'); else; opts.family = 'db4'; end
    [ denData ~ ]  = wDenoise(data, mask, opts);
    if isfield(opts, 'figdir')
	saveas(gcf,[opts.figdir,'wavDenoise_data.fig']);
	else
	saveas(gcf,[pwd,'\','wavDenoise_data.fig']);
	end
else
    [voxel_ts, coordinates] = grabTimeseries(data, mask);
    figure; plot(mean(voxel_ts,1),'b'); hold on;
    if isfield(opts,'method'); else; opts.method = 'movmean'; end
	if isfield(opts,'winsize'); else; opts.winsize = round(0.05*size(data,4))'; end %5 percent window

	voxel_ts = smoothdata(voxel_ts,2,opts.method);
    plot(mean(voxel_ts,1),'r');
    denData = reshape(data, [opts.xdim*opts.ydim*opts.zdim opts.dyn]);
    denData(coordinates,:) = voxel_ts;
    denData = reshape(denData, size(data));
    if isfield(opts, 'figdir')
	saveas(gcf,[opts.figdir,'denoise_data.fig']);
    else
	saveas(gcf,[pwd,'\','denoise_data.fig']);
	end
	clear voxel_ts coordinates
end

end

