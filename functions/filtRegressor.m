%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht,
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes.

function [innuisance outnuisance] = filtRegressor(nuisance, corrProbe, opts)
warning('off');
global opts
if isfield(opts,'motioncorr'); else; opts.motioncorr = 0.3; end

test1 = nuisance(1,:); test2 = nuisance(:,1);
if length(test1) > length(test2); nuisance = nuisance'; end; clear test1 test2
if iscolumn(corrProbe)==0; corrProbe = corrProbe'; end

%remove nuisance correlating with probe
autoCorr = abs(corr(rescale(corrProbe),nuisance)); 
autoCorr(autoCorr < opts.motioncorr) = 0; autoCorr(autoCorr > 0) = 1; %removes anything with more than weak correlation
%remove highly correlated nuisance regressors to preserve signal response
index = ([1:1:size(nuisance,2)]).*autoCorr; 
keep = index; keep(keep == 0) = [];
leave = index; leave = find(index == 0)
innuisance = nuisance; outnuisance = nuisance;

innuisance(:,keep) = [];
outnuisance(:,leave) = [];


end