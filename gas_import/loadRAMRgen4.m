% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <loadRAMRgen4: imports GEN4 RespirAct data >
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

function [corrvec_CO2,corrvec_O2, corrvec_RR ] = loadRAMRgen4(opts)
% This is the main function (used to call relevant sub-functions) to load
% the GEN4 respiratory data. This function can also perform outlier
% removal.
warning('off')
global opts;

if isfield(opts,'seqpath')
    cd(opts.seqpath)
else
    error('please specify location of RespirAct files in opts.seqpath')   
end

if isfield(opts,'extra'); else; opts.extra = 30; end                    %defines number of extra samples taken before and after the start and end of the events
if isfield(opts,'remOUT'); else; opts.remOUT = 0; end                   %remove outliers from breathing traces
if isfield(opts,'remOUTbh'); else; opts.remOUTbh = 0; end               %for breath hold, only outliers below the baseline are removed
if isfield(opts,'remOUTmethod'); else; opts.remOUTmethod = 'median'; end %'median' or 'mean' or 'quartiles' or 'grubbs' or 'gesd'


%Import EndTidal
file = ls('EndTidal.*');
filename = fullfile(opts.seqpath,file);
try
[MRTimes,DesiredPO2mmHg,DesiredPCO2mmHg,AchievablePO2mmHg,AchievablePCO2mmHg,PO2mmHg,PCO2mmHg,RestingPO2mmHg,RestingPCO2mmHg,PBarommHg,Inspiretimeseconds,Expiretimeseconds,Breathidx,TidalvolumemL,RespirationrateBPM,StartInspiresec,O2AdjustmentmmHg,CO2AdjustmentmmHg,G1TargetvolmL,G1FCO2,G1FO2,G2FCO2,G2FO2] = import_EndTidal(filename);
catch
    error('Check to make sure correct input folder is selected. fomr RAMR GEN4, this should be the folder that contains the participant and study sub-folders')
end

%import RGM
file = ls('RGM.*');
filename = fullfile(opts.seqpath,file);
[MRTimes1,PO2mmHg1,PCO2mmHg1,PBarommHg1,PMouthmmH2O,FlowMouthmLmin,FlowS1mLmin,FlowS2mLmin,BreathPhase] = import_RGM(filename);

%calculate picker derivative
brdiff = diff(MRTimes)


figure;
subplot(2,1,1)
plot(MRTimes1,PO2mmHg1)
hold on; scatter(MRTimes, PO2mmHg)
subplot(2,1,2)
plot(MRTimes, RespirationrateBPM)

%import events
file = ls('Events.*');
filename = fullfile(opts.seqpath,file);
try
[MRTimes3, CtrlRoomTimes, Event] = import_Events(filename);
catch
   error('Possible mismatch; rename events (i.e. prep, start session) to default values or adapt matlab code to match with your chosen event names') 
end
%import physiological parameters
file = ls('PhysioParams.*');
filename = fullfile(opts.seqpath,file);

[MRTimes2, ID, FRCmL, VdmL, TissuestoreO2mL, TissuestoreCO2mL, VO2mLmin, VCO2mLmin, QmLmin, hBconcentrationgdLBlood, Restingmetabolicscalefactor, ResponseReason] = import_Physio(filename);

% resample and realign the breathing trace and the Endtidal trace to have the same start and end and same sampling rate

[nxi,corrvec_CO2,corrvec_O2,nxi1,rawCO2,rawO2,corrvec_RR] = resampleEndtidalBreathing(MRTimes,PCO2mmHg,PO2mmHg,MRTimes1,PCO2mmHg1,PO2mmHg1,MRTimes3,RespirationrateBPM,Event,opts);


%% calculate RVT
end

