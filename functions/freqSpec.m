
function[ fSpec, meanSpec, freq ] = freqSpec(data, mask, opts)        
%this function returns the frequency spectrum series generated from 
%timeseries input data for a region defined by 'mask'. Additional return 
%variables are the average spectrum inthe ROI as well as the corresponding
%frequency values 

        if isfield(opts,'powerspec'); else; opts.powerspec = 0; end
        Fs = 1/opts.TR;
        %form frequency graphs
        [voxels,coordinates] = grabTimeseries(data, mask);
        [x,y,z,N] = size(data);
        freq = 0:Fs/N:Fs/2;
        r_xdft = fft(voxels,[],2);           % FFT
        r_xdft = r_xdft(:,1:N/2+1);           % take half of data
        if opts.powerspec
        r_psdx = (1/(Fs*N))*abs(r_xdft).^2;% power
        else
        r_psdx = (1/(Fs*N))*abs(r_xdft);% amplitude
        end
        r_psdx(2:end-1) = 2*r_psdx(2:end-1);  
        figure; set(gcf,'color','w');
        plot(freq,log(mean(r_psdx,1)));
        if opts.powerspec
        title('log of power spectrum'); 
        else
        title('log of amplitude spectrum');
        end
        hold on; 
        xline(0.01,'k--'); 
        xline(0.027,'k--');
        xline(0.073,'k--');
        xline(0.17,'k--');
        xline(0.23,'k--'); 
        hold off
        meanSpec = mean(r_psdx,1);
        fSpec = zeros([x*y*z,size(r_psdx,2)]);
        fSpec(coordinates,:) = r_psdx;
        fSpec = reshape(fSpec,[x,y,z,size(r_psdx,2)]); 
        disp('returning spectrum of mask region (NOT logtransformed)')
end