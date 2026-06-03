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
% This function resamples the end tidal values for gen4 respiract systems.
% extra defines the amount of extra time points before start and end of the
% sequence; there is also an option to remove outliers for breath hold data
% and normal data.
%
% Outputs:
%   xi    = resampled time vector for end-tidal traces
%   yi    = resampled end-tidal CO2
%   yyi   = resampled end-tidal O2
%   xi1   = resampled time vector for RGM traces
%   yi1   = resampled measured RGM CO2
%   yyi1  = resampled measured RGM O2
%   rri   = resampled respiration rate

%% Ensure column vectors

MRTimes    = MRTimes(:);
PCO2mmHg   = PCO2mmHg(:);
PO2mmHg    = PO2mmHg(:);

MRTimes1   = MRTimes1(:);
PCO2mmHg1  = PCO2mmHg1(:);
PO2mmHg1   = PO2mmHg1(:);

RR         = RR(:);

%% Set defaults

if isfield(opts,'figdir')
else
    if ispc
        opts.figdir = [pwd,'\']; 
    else
        opts.figdir = [pwd,'/'];
    end
end

%% Make sure MRTimes1 and MRTimes have the same start and end-point

% Use first and last value of MRTimes1 to put at beginning and end of MRTimes
startval = MRTimes1(1,1);
endval = MRTimes1(end,1);
nMRTimes = [startval; MRTimes; endval];

% Use first and last value of PCO2mmHg trace to extend it with
startval = PCO2mmHg(1,1);
endval = PCO2mmHg(end,1);
nPCO2mmHg = [startval; PCO2mmHg; endval];

% Use first and last value of PO2mmHg trace to extend it with
startval = PO2mmHg(1,1);
endval = PO2mmHg(end,1);
nPO2mmHg = [startval; PO2mmHg; endval];

% Use first and last value of RR trace to extend it with
startval = RR(1,1);
endval = RR(end,1);
nRR = [startval; RR; endval];

%% Resample entire PCO2mmHg and PO2mmHg trace to the TR of the acquisition

% Detect outlying points and replace with interpolated values
if opts.remOUT == 1

    TF = isoutlier(nPCO2mmHg,opts.remOUTmethod);

    figure(2); clf;
    set(gcf,'Color','w','Name','Outlier removal for breathhold');
    set(gcf,'Units','normalized','Position',[0.05 0.08 0.90 0.82]);

    subplot(2,1,1); 
    plot(nMRTimes,nPCO2mmHg,'.','Color',[0.00 0.85 1.00]);
    hold on; 
    plot(nMRTimes(TF),nPCO2mmHg(TF),'o', ...
        'Color',[1.00 0.00 0.85], ...
        'MarkerFaceColor',[1.00 0.00 0.85]);
    title('Outlier removal for CO2');
    ylabel('PCO2 [mmHg]');
    grid on;

    subplot(2,1,2); 
    plot(nMRTimes,nPO2mmHg,'.','Color',[1.00 0.00 0.85]);
    hold on; 
    plot(nMRTimes(TF),nPO2mmHg(TF),'o', ...
        'Color',[0.00 0.85 1.00], ...
        'MarkerFaceColor',[0.00 0.85 1.00]);
    title('Outlier removal for O2');
    ylabel('PO2 [mmHg]');
    xlabel('Time [s]');
    grid on;

    % Remove outliers lower than baseline mean
    ref_mean = mean(nPCO2mmHg(1:round(0.10*length(nPCO2mmHg))),1);
    idx = find(TF == 1);
    values = nPCO2mmHg(idx);

    if opts.remOUTbh
        values(values < ref_mean) = 0;
        filter = idx(values == 0);
    else
        filter = idx;
    end

    nPCO2mmHg(filter) = NaN; 
    nPO2mmHg(idx) = NaN;

    nPCO2mmHg = inpaint_nans(nPCO2mmHg); 
    nPO2mmHg = inpaint_nans(nPO2mmHg);

    figure(2);
    subplot(2,1,1); 
    plot(nMRTimes,nPCO2mmHg,'Color',[0.50 0.00 1.00],'LineWidth',1.2);

    subplot(2,1,2); 
    plot(nMRTimes,nPO2mmHg,'Color',[0.50 0.00 1.00],'LineWidth',1.2);

    saveas(gcf,[opts.figdir,'outlier_removal.fig']);

end

%% End-tidal values

[xi,yi] = resampletoTR(opts.TR,nMRTimes,nPCO2mmHg);   % yi  = resampled CO2
[xi,yyi] = resampletoTR(opts.TR,nMRTimes,nPO2mmHg);  % yyi = resampled O2
[~,rri] = resampletoTR(opts.TR,nMRTimes,nRR);        % rri = resampled respiration rate

%% RGM values

% Keep original measured RGM values for returned outputs.
% These are not used quantitatively in the lower-row diagnostic y-axis.
[xi1,yi1] = resampletoTR(opts.TR,MRTimes1,PCO2mmHg1); % yi1  = resampled measured RGM CO2
[~,yyi1] = resampletoTR(opts.TR,MRTimes1,PO2mmHg1);   % yyi1 = resampled measured RGM O2

%% Prepare scaled RGM breathing waveform for visualization only

% The raw RGM y-values should not be interpreted as end-tidal mmHg here.
% We use the RGM trace only to show the inhale/exhale waveform shape.
% Therefore, scale the RGM waveform into the y-range of the corresponding
% end-tidal values for plotting only.

PCO2_RGM_scaled = scaleTraceToReferenceRange(PCO2mmHg1, PCO2mmHg);
PO2_RGM_scaled  = scaleTraceToReferenceRange(PO2mmHg1,  PO2mmHg);

%% Visual check

% Colors: cyan / magenta / purple style
cyan          = [0.00 0.85 1.00];
cyan_dark     = [0.00 0.45 0.60];
cyan_light    = [0.65 0.95 1.00];

magenta       = [1.00 0.00 0.85];
magenta_dark  = [0.65 0.00 0.55];
magenta_light = [1.00 0.70 0.95];

purple        = [0.45 0.00 1.00];
purple_dark   = [0.25 0.00 0.65];

yellow        = [1.00 0.75 0.00];
grey          = [0.35 0.35 0.35];

figure(1); clf;
set(gcf,'Color','w','Name','RespirAct GEN4 resampling diagnostic');
set(gcf,'Units','normalized','Position',[0.02 0.05 0.96 0.86]);

% -------------------------------------------------------------------------
% Row 1: End-tidal values
% -------------------------------------------------------------------------

subplot(3,3,1)
plot(nMRTimes,nPCO2mmHg,'-','Color',cyan_dark,'LineWidth',1.2);
hold on;
plot(MRTimes,PCO2mmHg,'o', ...
    'Color',cyan_dark, ...
    'MarkerFaceColor','w', ...
    'MarkerSize',4);
title('End-tidal CO2 values');
ylabel('PCO2 [mmHg]');
xlabel('End-tidal time [s]');
legend({'Extended ET CO2','Original ET CO2'},'Location','best');
grid on;

subplot(3,3,2)
plot(nMRTimes,nPO2mmHg,'-','Color',magenta_dark,'LineWidth',1.2);
hold on;
plot(MRTimes,PO2mmHg,'o', ...
    'Color',magenta_dark, ...
    'MarkerFaceColor','w', ...
    'MarkerSize',4);
title('End-tidal O2 values');
ylabel('PO2 [mmHg]');
xlabel('End-tidal time [s]');
legend({'Extended ET O2','Original ET O2'},'Location','best');
grid on;

subplot(3,3,3)
yyaxis left
plot(nMRTimes,nPCO2mmHg,'-','Color',cyan_dark,'LineWidth',1.2);
ylabel('PCO2 [mmHg]');

yyaxis right
plot(nMRTimes,nPO2mmHg,'-','Color',magenta_dark,'LineWidth',1.2);
ylabel('PO2 [mmHg]');

title('End-tidal CO2 and O2');
xlabel('End-tidal time [s]');
grid on;

% -------------------------------------------------------------------------
% Row 2: Resampled end-tidal values
% -------------------------------------------------------------------------

subplot(3,3,4)
plot(xi,yi,'-','Color',cyan_dark,'LineWidth',1.4);
hold on;
plot(xi,yi,'.','Color',cyan_dark,'MarkerSize',8);
title('Resampled end-tidal CO2');
ylabel('PCO2 [mmHg]');
xlabel('Resampled time [s]');
grid on;

subplot(3,3,5)
plot(xi,yyi,'-','Color',magenta_dark,'LineWidth',1.4);
hold on;
plot(xi,yyi,'.','Color',magenta_dark,'MarkerSize',8);
title('Resampled end-tidal O2');
ylabel('PO2 [mmHg]');
xlabel('Resampled time [s]');
grid on;

subplot(3,3,6)
yyaxis left
plot(xi,yi,'-','Color',cyan_dark,'LineWidth',1.4);
ylabel('PCO2 [mmHg]');

yyaxis right
plot(xi,yyi,'-','Color',magenta_dark,'LineWidth',1.4);
ylabel('PO2 [mmHg]');

title('Resampled ET CO2 and O2');
xlabel('Resampled time [s]');
grid on;

% -------------------------------------------------------------------------
% Row 3: RGM breathing waveform shape with raw ET values overlaid
% -------------------------------------------------------------------------
% Important:
% The RGM y-values are NOT used quantitatively here.
% The waveform is scaled into the ET y-range only to show breathing phase.
% The yellow points are the actual measured end-tidal values.

subplot(3,3,7)
plot(MRTimes1,PCO2_RGM_scaled,'-','Color',cyan_light,'LineWidth',0.8);
hold on;
plot(MRTimes1,PCO2_RGM_scaled,'.','Color',cyan,'MarkerSize',3);

plot(MRTimes,PCO2mmHg,'o', ...
    'Color',purple_dark, ...
    'MarkerFaceColor',yellow, ...
    'MarkerSize',5);

title('RGM waveform shape + raw end-tidal CO2');
ylabel('PCO2 [mmHg]');
xlabel('Time [s]');
legend({'RGM waveform scaled to ET range','RGM samples','Raw ET CO2'}, ...
    'Location','best');
grid on;

subplot(3,3,8)
plot(MRTimes1,PO2_RGM_scaled,'-','Color',magenta_light,'LineWidth',0.8);
hold on;
plot(MRTimes1,PO2_RGM_scaled,'.','Color',magenta,'MarkerSize',3);

plot(MRTimes,PO2mmHg,'o', ...
    'Color',purple_dark, ...
    'MarkerFaceColor',yellow, ...
    'MarkerSize',5);

title('RGM waveform shape + raw end-tidal O2');
ylabel('PO2 [mmHg]');
xlabel('Time [s]');
legend({'RGM waveform scaled to ET range','RGM samples','Raw ET O2'}, ...
    'Location','best');
grid on;

subplot(3,3,9)
plot(MRTimes1,PCO2_RGM_scaled,'-','Color',cyan,'LineWidth',1.0);
hold on;
plot(MRTimes1,PO2_RGM_scaled,'-','Color',magenta,'LineWidth',1.0);
title('Scaled RGM waveform shapes');
ylabel('Scaled to ET ranges');
xlabel('Time [s]');
legend({'CO2 RGM waveform shape','O2 RGM waveform shape'},'Location','best');
grid on;

%% Save diagnostic figure

try
    saveas(gcf,[opts.figdir,'respiract_resampling_diagnostic.fig']);
catch
    warning('Could not save RespirAct diagnostic figure.');
end

end

%% Local helper function

function yScaled = scaleTraceToReferenceRange(yRaw, yRef)

% Scale yRaw into the y-range of yRef for visualization only.
% This preserves waveform shape but discards the original RGM y-units.

yRaw = yRaw(:);
yRef = yRef(:);

validRaw = isfinite(yRaw);
validRef = isfinite(yRef);

yScaled = yRaw;

if sum(validRaw) < 2 || sum(validRef) < 2
    return;
end

rawMin = prctile(yRaw(validRaw),1);
rawMax = prctile(yRaw(validRaw),99);

refMin = min(yRef(validRef));
refMax = max(yRef(validRef));

% Add small margins so waveform does not sit exactly on the plot borders
refRange = refMax - refMin;
refMin = refMin - 0.05 * refRange;
refMax = refMax + 0.05 * refRange;

if rawMax == rawMin || refMax == refMin
    yScaled(validRaw) = mean([refMin refMax]);
else
    yScaled(validRaw) = (yRaw(validRaw) - rawMin) ./ (rawMax - rawMin);
    yScaled(validRaw) = yScaled(validRaw) .* (refMax - refMin) + refMin;
end

end