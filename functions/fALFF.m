function [ALFF_map fALFF_map zALFF_map zfALFF_map] = fALFF(data, mask, refmask, opts)                  
%Script written by ALEX BHOGAL a.bhogal@umcutrecht.nl
%   Detailed explanation goes here
%timeseries data is used to generate ALFF and fALFF maps based on specified
%frequency band. 
%see: An improved approach to detection of amplitude of low-frequency fluctuation (ALFF) for resting-state fMRI: Fractional ALFF
%doi: 10.1016/j.jneumeth.2008.04.012
global opts
opts.ALFFdir = [opts.resultsdir,'ALFF/']; mkdir(opts.ALFFdir);
[xx yy zz N] = size(data);
[refdata] = meanTimeseries(data, mask);

opts.voxelsize = opts.headers.map.dime.pixdim(2:4)
Fs = 1/opts.TR;
dF = Fs/N;
xdft = fft(refdata);
xdft = xdft(1:N/2+1);
%reference signal
psdx = abs(xdft)/N; %square this to get the power (ampl^2 = pwer)
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
psdx = abs(xdft)/N; %square this to get the power (ampl^2 = pwer)
psdx(2:end-1) = 2*psdx(2:end-1);

vALFF = sum(sqrt(psdx(:,Lowidx:Highidx)),2); %ALFF in freq band of interest
fvALFF = sum(sqrt(psdx),2); %whole spectrum reference ALFF

%generate ALFF map and z-transformed ALFF map
ALFF_map = zeros(size(mask)); ALFF_map = ALFF_map(:); zALFF_map = ALFF_map;
ALFF_map(coordinates,1) =  vALFF/rALFF; %divide by global mean ALFF - see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3902859/
zALFF_map(coordinates,1) = (vALFF - mean(vALFF))/std(vALFF); %calculate z-score
ALFF_map = reshape(ALFF_map,size(mask)); zALFF_map = reshape(zALFF_map,size(mask));
ALFF_map = smthData( ALFF_map, double(mask), opts); zALFF_map = smthData(zALFF_map, double(mask), opts);
limits = string(opts.fpass);
%save ALFF
name = join(['ALFF_map_',limits(1),'_',limits(2),'.nii.gz']);
name(name == ' ') = [];
saveImageData(ALFF_map,opts.headers.map,opts.ALFFdir,char(name),64)
%save ALFF z-map
name = join(['zALFF_map_',limits(1),'_',limits(2),'.nii.gz']);
name(name == ' ') = [];
saveImageData(zALFF_map,opts.headers.map,opts.ALFFdir,char(name),64)

%generate ALFF map and z-transformed ALFF map
fALFF_map = zeros(size(mask)); fALFF_map = fALFF_map(:); zfALFF_map = fALFF_map;
fALFF_map(coordinates,1) =  vALFF./fvALFF;
zfALFF_map(coordinates,1) = (fALFF_map(coordinates,1) - mean(fALFF_map(coordinates,1)))/std(fALFF_map(coordinates,1));
fALFF_map = reshape(fALFF_map,size(mask)); zfALFF_map = reshape(zfALFF_map,size(mask));
fALFF_map = smthData( fALFF_map, double(mask), opts); zfALFF_map = smthData(zfALFF_map, double(mask), opts);
%save fALFF
name = join(['fALFF_map_',limits(1),'_',limits(2),'.nii.gz']);
name(name == ' ') = [];
saveImageData(fALFF_map,opts.headers.map,opts.ALFFdir,char(name),64)
%save fALFF z-map
name = join(['zfALFF_map_',limits(1),'_',limits(2),'.nii.gz']);
name(name == ' ') = [];
saveImageData(zfALFF_map,opts.headers.map,opts.ALFFdir,char(name),64)
end
