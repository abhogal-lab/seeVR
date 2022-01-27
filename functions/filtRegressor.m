% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <filtRegressor: correlates a series of regressors with a reference signal >
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
% This function correlates a series of nuisance regressors (typically 
% motion parameters derived from realignment) with the input probe. Based
% on the correlation threshold (opts.motioncorr, nuisance regressors are
% seperated.
%
% nuisance: an array of nuisance regressors (i.e. motion and derivatives,
% HR, respiratory etc.)
%
% probe: timeseries to which nuisance parameters are correlated
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.motioncorr
%
% innuisance: regressors having a correlation lower than the threshold
%
% outnuisance: regressors having a correlation higher than the threshold
function [innuisance outnuisance] = filtRegressor(nuisance, probe, opts)


warning('off');
global opts
if isfield(opts,'motioncorr'); else; opts.motioncorr = 0.3; end

test1 = nuisance(1,:); test2 = nuisance(:,1);
if length(test1) > length(test2); nuisance = nuisance'; end; clear test1 test2
if iscolumn(probe) == 0; probe = probe'; end

%remove nuisance correlating with probe
autoCorr = abs(corr(probe,nuisance)); 
autoCorr(autoCorr < opts.motioncorr) = 0; autoCorr(autoCorr > 0) = 1; %removes anything with more than weak correlation
%remove highly correlated nuisance regressors to preserve signal response
index = ([1:1:size(nuisance,2)]).*autoCorr; 
keep = index; keep(keep == 0) = [];
leave = index; leave = find(index == 0);

innuisance = nuisance; outnuisance = nuisance;

innuisance(:,keep) = [];
if isempty(innuisance)
    disp('Motion threshold is set too low. Increase opts.motioncorr and try again')
    return
end
outnuisance(:,leave) = [];

% Plot
if isempty(outnuisance)
cols = 3; else; cols = 4; 
end
xdata = opts.TR:opts.TR:length(probe)*opts.TR;
figure; 
subplot(cols,1,1); plot(xdata, probe); 
title('reference signal'); xlabel('time(s)'); 
xlim([0 xdata(end)]);
set(gcf, 'Units', 'pixels', 'Position', [200, 500, 600, 800]);
subplot(cols,1,2); plot(xdata, nuisance); 
title('motion parameters & derivatives'); xlabel('time(s)'); 
xlim([0 xdata(end)]);
subplot(cols,1,3); plot(xdata, innuisance); 
title('nuisance < r-threshold'); xlabel('time(s)');
xlim([0 xdata(end)]);
if ~isempty(outnuisance)
subplot(cols,1,4); plot(xdata, outnuisance); 
title('nuisance > r-threshold'); xlabel('time(s)');
xlim([0 xdata(end)]);
end
end