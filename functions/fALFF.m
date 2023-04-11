% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <fALFF: calculates (fractional) amplitude of low frequency fluctuations ((f)ALFF) >
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
% Timeseries data is used to generate ALFF and fALFF maps based on specified
% frequency bands. For details see:
% An improved approach to detection of amplitude of low-frequency fluctuation (ALFF) for resting-state fMRI: Fractional ALFF
%doi: 10.1016/j.jneumeth.2008.04.012
%
% data: input timeseries data (i.e. 4D BOLD MRI dataset)
%
% mask: binary mask defining voxels of interest
%
% refmask: binary mask that defines the reference signals of interest
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.fpass, opts.resultsdir
%
% ALFF_map: ALFF map defined by opts.fpass
%
% fALFF_map: fALFF map defined by opts.fpass
%
% zALFF_map: z-transformed ALFF map defined by the standard deviation of
% ALFF values defined by the input mask
%
% zfALFF_map: z-transformed fALFF map defined by the standard deviation of
% fALFF values defined by the input mask
function [ALFF_map fALFF_map zALFF_map zfALFF_map] = fALFF(data, mask, refmask, opts)
global opts

if isfield(opts,'niiwrite'); else; opts.niiwrite = 0; end                  %depending on how data is loaded this can be set to 1 to use native load/save functions

opts.ALFFdir = fullfile(opts.resultsdir,'ALFF'); mkdir(opts.ALFFdir);

[xx yy zz N] = size(data);
[refdata] = meanTimeseries(data, mask);

Fs = 1/opts.TR;
dF = Fs/N;
xdft = fft(refdata);
xdft = xdft(1:N/2+1);
%reference signal
psdx = (abs(xdft).*2)/N; %square this to get the power (ampl^2 = pwer)
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(data):Fs/2;

Lowf = opts.fpass(1); Highf = opts.fpass(2);

%find idx for low & high freq
[Lowval,Lowidx]=min(abs(freq-Lowf));
[Highval,Highidx]=min(abs(freq-Highf));

rALFF = sum(sqrt(psdx(1,Lowidx:Highidx))); %whole brain reference ALFF see: https://pubmed.ncbi.nlm.nih.gov/16919409/

%voxel-wise ALFF
[voxels coordinates] = grabTimeseries(data, mask);
xdft = fft(voxels,[],2);
xdft = xdft(:,1:N/2+1);
psdx = (abs(xdft).^2)/N; %square this to get the power (ampl^2 = pwer)
psdx(2:end-1) = 2*psdx(2:end-1); %not sure yet is factor 2 is needed

vALFF = sum(sqrt(psdx(:,Lowidx:Highidx)),2); %ALFF in freq band of interest
fvALFF = sum(sqrt(psdx),2); %whole spectrum reference ALFF

%generate ALFF map and z-transformed ALFF map
ALFF_map = zeros(size(mask)); ALFF_map = ALFF_map(:); zALFF_map = ALFF_map;
ALFF_map(coordinates,1) =  vALFF/rALFF; %divide by global mean ALFF - see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3902859/
zALFF_map(coordinates,1) = (vALFF - mean(vALFF))/std(vALFF); %calculate z-score
ALFF_map = reshape(ALFF_map,size(mask)); zALFF_map = reshape(zALFF_map,size(mask));

limits = string(opts.fpass);

%save ALFF
delimeter = {'_','_','_'};
name = join(['ALFF_map',limits(1),limits(2),'.nii.gz'],delimeter);
if opts.niiwrite
    niftiwrite(cast(ALFF_map,opts.mapDatatype),[opts.ALFFdir,name],opts.info.map);
else
    saveImageData(ALFF_map,opts.headers.map,opts.ALFFdir,char(name),64)
end

%save ALFF z-map
name = join(['zALFF_map',limits(1),limits(2),'.nii.gz'],delimeter);
if opts.niiwrite
    cd(opts.ALFFdir)
    niftiwrite(cast(zALFF_map, opts.mapDatatype),name,opts.info.map);
else
    saveImageData(zALFF_map,opts.headers.map,opts.ALFFdir,char(name),64)
end

%generate ALFF map and z-transformed ALFF map
fALFF_map = zeros(size(mask)); fALFF_map = fALFF_map(:); zfALFF_map = fALFF_map;
fALFF_map(coordinates,1) =  vALFF./fvALFF;
zfALFF_map(coordinates,1) = (fALFF_map(coordinates,1) - mean(fALFF_map(coordinates,1)))/std(fALFF_map(coordinates,1));
fALFF_map = reshape(fALFF_map,size(mask)); zfALFF_map = reshape(zfALFF_map,size(mask));

%save fALFF
name = join(['fALFF_map',limits(1),limits(2),'.nii.gz'],delimeter);
if opts.niiwrite
    niftiwrite(cast(fALFF_map, opts.mapDatatype),name,opts.info.map);
else
    saveImageData(fALFF_map,opts.headers.map,opts.ALFFdir,char(name),64)
end

%save fALFF z-map
name = join(['zfALFF_map',limits(1),limits(2),'.nii.gz'],delimeter);
if opts.niiwrite
    niftiwrite(cast(zfALFF_map, opts.mapDatatype),name,opts.info.map);
else
    saveImageData(zfALFF_map,opts.headers.map,opts.ALFFdir,char(name),64)
end
end
