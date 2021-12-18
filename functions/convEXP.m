% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <convHRF: convolved an input signal with a double-gamma hemodynamic response function (HRF) >
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
% ************************************************************************

function [HRF, HRF_probe] = convEXP(probe,opts)
warning('off')
global opts;
%set defaults
if isfield(opts,'verbose'); else; opts.verbose = 0; end %turn plotting on or off
if isfield(opts,'disp'); else; opts.disp = (1:.5:30); end %dispersion parameter - set a default vector
if isfield(opts,'pad'); else; opts.pad = 1; end %pad before fft
if isfield(opts,'padfront'); else; opts.padfront = 0; end %pad before fft

input_probe = rescale(probe);

if opts.pad
    pad = 2^nextpow2(length(probe));
else
    pad = 0;
end
if opts.padfront && opts.pad
    probe = zeros(length(input_probe)+pad,1);
    probe(1:pad) = 0;% input_probe(1,1);
    probe(pad+1:end,1) = input_probe;
else
    probe = zeros(length(input_probe)+2*pad,1);
    probe(1:pad) = 0; %input_probe(1,1);
    probe(pad+1:end-pad,1) = input_probe;
    probe(end-pad+1:end,1) = 0; %input_probe(1,end);
end

% generate HRF function
xdata = (1:1:size(probe,1));
t = (0:1:length(xdata)-1);
HRF_probe = zeros([length(opts.disp)+1 size(probe,1)]);
HRF_probe(1,:) = probe/trapz(probe);

for ii=1:length(opts.disp)   
HRF(ii,:) = exp(-xdata/(opts.disp(ii)));
HRF_freq = fft(HRF(ii,:)'); %HRF_freq/trapz(HRF_freq);
newprobe_freq = fft(probe); %newprobe_freq/trapz(newprobe_freq);
tmp = abs(ifft(HRF_freq.*newprobe_freq));
tmp = tmp/trapz(tmp); %normalize HRF
HRF_probe(ii+1,:) = tmp;%-min(tmp);
end

if opts.padfront && opts.pad
    HRF_probe(:,1:pad) = [];
else
    HRF_probe(:,1:pad) = []; HRF_probe(:,end-pad+1:end) = [];  %remove padding
end
for ii=1:size(HRF_probe,1); HRF_probe(ii,:) = rescale(HRF_probe(ii,:)); end
% remove linearly independent components
[Xsub,idx]=licols(HRF_probe');
XHRFidx = HRFidx(idx,:);
HRF_probe = Xsub';
HRFidx = XHRFidx;

%remove artificial trend
if opts.detrendHRF
    %use first and last 10% of data pts (improve this later)
    HRF_probe = detrendHRF(HRF_probe, [50 size(HRF_probe,2)-20]);
end

for ii=1:size(HRF_probe,1); HRF_probe(ii,:) = rescale(HRF_probe(ii,:)); end
if opts.verbose
    figure;
    %customMap = viridis(size(HRF,1)); customMap = flip(customMap,1);
    customMap = colormap(flip(brewermap(size(HRF,1),'Spectral')));
    for ii=1:size(HRF,1)
        subplot(1,2,1); plot(t,HRF(ii,:)', 'Color', customMap(ii,:)); ylim([-0.2 1]); xlim([0 200]); title('exponential HRF'); hold on;
        subplot(1,2,2); plot(HRF_probe(ii+1,:)', 'Color', customMap(ii,:)); xlim([0 size(HRF_probe,2)]); hold on;
    end
    plot( HRF_probe(1,:), 'k', 'LineWidth', 2);  title('HRF probe');
    set(gcf, 'Units', 'pixels', 'Position', [200, 500, 600, 160]);
    if isfield(opts,'figdir')
        saveas(gcf,[opts.figdir,'EXP_fnc.fig']);
    else
        if ispc
            saveas(gcf,[pwd,'\','EXP_fnc.fig']);
        else
            saveas(gcf,[pwd,'/','EXP_fnc.fig']);
        end
    end
end
end