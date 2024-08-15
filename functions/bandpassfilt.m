% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <bandpassfilt: performs a bandpass filter on input data >
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
% data: input timeseries data (i.e. 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% refmask: binary mask that defines the reference signals of interest
%
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.fpass opts.headers.map, opts.info.map, opts,niiwrite, opts.resultsdir
%
%
% ROI defined by refmask
%
% bpData: a bandpassed version of the input data%
function [bpData] = bandpassfilt(data, mask, opts)
warning('off');
global opts;
t = cputime;

% data = normTimeseries(data, mask);
mean_image = mean(data, 4);
data = demeanData(data, mask);

mask = logical(mask);

if isfield(opts,'fpass'); else; opts.fpass = [0.000001 0.1164]; end  %default frequency band see doi: 10.1148/radiol.2021203568
    
[xx yy zz N] = size(data);
data = double(data);

%get voxels
[voxels_ts coordinates] = grabTimeseries(data, mask);

Fs = 1/opts.TR;
t = 0:1/Fs:1-1/Fs;
N = size(voxels_ts,2);
dF = Fs/N;
freq = 0:Fs/length(mean_ts):Fs/2;
Lowf = opts.fpass(1); Highf = opts.fpass(2); %Hz

%check whether signal processing toolbox is available
if license('test','signal_toolbox') == 1
    % if opts.filter_order is specified
    if isfield(opts,'filter_order')
        [b,a] = butter(opts.filter_order,2*[Lowf, Highf]/Fs);
        BP_ref = filtfilt(b,a,(mean_ts));
    else
        % else find optimal filter (automated)
        SSres = []; BP = [];
        for forder = 1:4
            [b,a] = butter(forder,2*[Lowf, Highf]/Fs);
            BP(forder,:) = filtfilt(b,a,(mean_ts));
            SSres(forder,:)= sum((mean_ts - BP(forder,:).^2),2);
        end
        %select closest filter
        [~,I] = min(SSres);
        [b,a] = butter(I,2*[Lowf, Highf]/Fs);
        BP_ref = filtfilt(b,a,(mean_ts));
        opts.filter_order = I;
    end
    
    %bandpass signals of interest contained within mask
    BP_V = zeros(size(voxels_ts));
    parfor ii=1:size(BP_V,1)
        BP_V(ii,:) = filtfilt(b,a,voxels_ts(ii,:));
    end
    
else
    isplot = 0;

    %bandpass signals of interest contained within mask
    BP_V = zeros(size(voxels_ts));
    BP_V = bpfilt(voxels_ts,Lowf, Highf, isplot);
end

%regress whole brain & reference signals against band passed reference

bpData = zeros(size(data));
bpData = reshape(bpData,[xx*yy*zz N]);
bpData(coordinates, :) = BP_V;
bpData = reshape(bpData,size(data)) + mean_image;

disp(['finished bandpass filtering of input data in: ',int2str(cputime-t),' seconds']);

end



