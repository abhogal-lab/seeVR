%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht,
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes.

function [denData] = denoiseData(data,mask,opts)

if license('test', 'wavelet_toolbox')
    if isfield(opts,'wdlevel'); else; opts.wdlevel = 3; end
    if isfield(opts,'family'); else; opts.family = 'db4'; end
    [ denData ~ ]  = wDenoise(data, mask, opts);
    if isfield(opts, 'figdir')
        saveas(gcf,[opts.figdir,'wavDenoise_data.fig']);
    else
        if ispc
            saveas(gcf,[pwd,'\','wavDenoise_data.fig']);
        else
            saveas(gcf,[pwd,'/','wavDenoise_data.fig']);
        end
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
        if ispc
        saveas(gcf,[pwd,'\','denoise_data.fig']);
        else
        saveas(gcf,[pwd,'/','denoise_data.fig']);
        end
        end
    clear voxel_ts coordinates
end

end

