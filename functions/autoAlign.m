% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <autoAlign: uses a 2-step GLM to search for optimal alignment between signals >
% This function is inspired by the work presented by Liu et al., (2022) https://doi.org/10.1371/journal.pone.0274220
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
% probe1: input to which fMRI time-series is to be aligned with (usually
% end-tidal CO2 trace
%
% probe2: optional second signal for alignment (usually end-tidal O2 or
% some other concurrently acquired physiological signal)
%
% TS: input time-series for alignment (usually grey matter or whole brain
% fMRI signal
%
% *note all signals should be interpolated to the same timing - i.e. to the
% TR of the fMRI scan. The functional ; resampletoTR.m can be used for this

function [probe1a,probe2a] = autoAlign(probe1,probe2, TS)

varCheck = nargin;
switch varCheck
    case 3
        if iscolumn(probe1); else probe1 = probe1'; end
        if iscolumn(probe2); else probe2 = probe2'; end
        if iscolumn(TS); else TS = TS'; end
        TS = rescale(TS);
    case 2
        if iscolumn(probe1); else probe1 = probe1'; end
        TS = probe2;
        if iscolumn(TS); else TS = TS'; end
        TS = rescale(TS);
end

% generate sliding window
regr = [];
for ii=1:length(probe1)-length(TS)
    regr(:,ii) = probe1(ii:ii+length(TS)-1,1);
end


%perform initial regression
regr_coef = zeros([size(regr,2) 2]);
for ii=1:size(regr,2)
    A = regr(:,ii);
    C = [ones([length(A) 1]) A];
    regr_coef(ii,:)= C\TS;
end

[M,I] = max(regr_coef(:,2));

try
    
    newprobe1 = probe1(I-8:I+length(TS)+7,1);
    %interpolate trace and data by a factor 10
    newprobe1 = interp(newprobe1,10);
    TS1 = interp(TS,10);
    
    if nargin == 3
        newprobe2 = probe2(I-8:I+length(TS)+7,1);
        newprobe2 = interp(newprobe2,10);
    else
        probe2a = []; %return empty variable
    end
    
    %setup new regr using sliding window
    regr = [];
    for ii=1:length(newprobe1)-length(TS1)
        regr(:,ii) = newprobe1(ii:ii+length(TS1)-1,1);
    end
    
    regr_coef = zeros([size(regr,2) 2]);
    for ii=1:size(regr,2)
        A = regr(:,ii);
        C = [ones([length(A) 1]) A];
        regr_coef(ii,:)= C\TS1;
    end
    
    %calculate residuals
    
    SSE = zeros([1 size(regr_coef,1)]); SST = zeros([1 size(regr_coef,1)]);  T = zeros([1 size(regr_coef,1)]);
    
    %predictions
    Y = regr_coef(:,2)*TS1' + regr_coef(:,2);
    %data
    X = repmat(TS1,1,size(regr_coef,1));
    
    for ii=1:size(X,2)
        SSE(1,ii) = (norm(X(:,ii) - Y(ii,:)'))^2;
        SST(1,ii) = (norm(X(:,ii)-mean(X(:,ii))))^2;
    end
    
    R2 = 1 - SSE./SST;
    x = [1:1:length(SSE)];
    p = polyfit(x,SSE,4);
    
    % Evaluate the fitted polynomial p and plot:
    f = polyval(p,x);
    figure; subplot(2,1,1)
    plot(x,SSE,'o',x,f,'-')
    legend('SSE','4th order polynomial fit');
    title('Refined shift search')
    
    [M,I] = min(f);
    if nargin == 3
        probe1a = newprobe1(I:I+length(TS1)-1,1);
        probe2a = newprobe2(I:I+length(TS1)-1,1);
        %downsample
        probe1a = downsample(probe1a, 10);
        probe2a = downsample(probe2a, 10);
        TS1 =downsample(TS1, 10);
        subplot(2,1,2); plot(rescale(probe1a),'m'); hold on; plot(rescale(TS1),'k'); hold on; plot(rescale(probe2a),'c')
    else
        probe1a = newprobe1(I:I+length(TS1)-1,1);
        %downsample
        probe1a = downsample(probe1a, 10);
        TS1 = downsample(TS1, 10);
        subplot(2,1,2); plot(rescale(probe1a),'m'); hold on; plot(rescale(TS1),'k')
    end
catch
    error('check input data or use manual alignment')
end

end

