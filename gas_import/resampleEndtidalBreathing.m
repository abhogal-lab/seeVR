% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <resampleEndtidalBreathing: resampled breathing data to the time-step defined by TR >
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

function [xi,yi,yyi,xi1,yi1,yyi1,rri] = resampleEndtidalBreathing(MRTimes,PCO2mmHg,PO2mmHg,MRTimes1,PCO2mmHg1,PO2mmHg1,MRTimes3,RR,Event,opts)
% This function resamples the end tidal values for gen4 respiract systems
% extra defines the amount of extra time points before start and end of the
% sequence; there is also an option to remove outliers for breath hold data
% and normal data.

%% make sure MRTimes1 and MRTimes have the same start and end-point
%use first and last value of MRTimes1 to put at begin and end of MRTimes
if isfield(opts,'figdir'); else
    if ispc
    opts.figdir = [pwd,'\']; 
else
    opts.figdir = [pwd,'/'];
    end
end

startval = MRTimes1(1,1);
endval = MRTimes1(end,1);
nMRTimes = [startval; MRTimes; endval];

%use first and last value of PCO2mmHg trace to extend it with
startval = PCO2mmHg(1,1);
endval = PCO2mmHg(end,1);
nPCO2mmHg = [startval; PCO2mmHg; endval];

%use first and last value of PO2mmHg trace to extend it with
startval = PO2mmHg(1,1);
endval = PO2mmHg(end,1);
nPO2mmHg = [startval; PO2mmHg; endval];

%use first and last value of PO2mmHg trace to extend it with
startval = RR(1,1);
endval = RR(end,1);
nRR = [startval; RR; endval];

figure(1); %visual check
subplot(2,3,1)
plot(nMRTimes,nPCO2mmHg,'b'); title('Endtidal CO2 values'); hold on; plot(MRTimes,PCO2mmHg,'r');
subplot(2,3,4) %visual check
plot(nMRTimes,nPO2mmHg,'b'); title('Endtidal O2 values'); hold on; plot(MRTimes,PO2mmHg,'r');

%% resample entire PCO2mmHg and PO2mmHg trace to the TR of the acquisition

% detect outlying points and replace with interpolated values
if opts.remOUT == 1
TF = isoutlier(nPCO2mmHg,opts.remOUTmethod);
figure(2); title('Outlier removal for breathold');
subplot(2,1,1); plot(nMRTimes,nPCO2mmHg,'.b');
hold on; plot(nMRTimes(TF),nPCO2mmHg(TF),'or');
subplot(2,1,2); plot(nMRTimes,nPO2mmHg,'.b');
hold on; plot(nMRTimes(TF),nPO2mmHg(TF),'or');
% remove ouliers lower than baseline mean
ref_mean = mean(nPCO2mmHg(1:round(0.10*length(nPCO2mmHg))),1);
idx = find(TF == 1);
values = nPCO2mmHg(idx);
if opts.remOUTbh
values(values < ref_mean) = 0;
filter = idx(values == 0);
else
filter = idx;
end
nPCO2mmHg(filter) = NaN; nPO2mmHg(idx) = NaN;
nPCO2mmHg = inpaint_nans(nPCO2mmHg); nPO2mmHg = inpaint_nans(nPO2mmHg); %fill removed outlier point
figure(2);
subplot(2,1,1); plot(nMRTimes,nPCO2mmHg,'m');
subplot(2,1,2); plot(nMRTimes,nPO2mmHg,'m');
saveas(gcf,[opts.figdir,'outlier_removal.fig']);
end

%%% End-tidal values %%%
[xi,yi] = resampletoTR(opts.TR,nMRTimes,nPCO2mmHg); %yi = resampled CO2
[xi,yyi] = resampletoTR(opts.TR,nMRTimes,nPO2mmHg); %yyi = resampled O2
[ri,rri] = resampletoTR(opts.TR,nMRTimes,nRR); %rri = resampled O2

%visual check
figure(1);
subplot(2,3,2); plot(nMRTimes,nPCO2mmHg); title('raw trace'); hold on; plot(nMRTimes,nPCO2mmHg,'.'); 
subplot(2,3,5); plot(xi,yi); title('resampled trace'); hold on; plot(xi,yi,'.');

%%% RGM values %%%
[xi1,yi1] = resampletoTR(opts.TR,MRTimes1,PCO2mmHg1); %yi = resampled CO2
[xil,yyi1] = resampletoTR(opts.TR,MRTimes1,PO2mmHg1); %yyi = resampled O2

%visual check
figure(1);
subplot(2,3,3); plot(MRTimes1,PCO2mmHg1); title('raw trace'); hold on; plot(MRTimes1,PCO2mmHg1,'.'); 
subplot(2,3,6); plot(xi1,yi1); title('resampled trace'); hold on; plot(xi1,yi1,'.');



end
