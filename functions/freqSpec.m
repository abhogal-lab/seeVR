% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <freqSpec: returns the frequency amplitude (or power) spectrum of input timeseries data >
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
% this function takes input data and uses a fourier analysis to generate
% the frequency spectrum information. If opts.powerspec = 1 is supplied in
% the opts structure, then the power spectrum will be returned.
%
% data: input timeseries data (i.e. 3D (slice) or 4D (volume) MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.powerspec
%
% fSpec: the frequency spectrum data derived from the input data. If
% opts.powerspec = 1, this variable contains the power spectrum.
%
% meanSpec: the average amplitude/power spectrum of the ROI defined by the
% mask
%
% freq: the vector of frequency values
function[ fSpec, meanSpec, freq ] = freqSpec(data, mask, opts)        

warning('off')
global opts
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