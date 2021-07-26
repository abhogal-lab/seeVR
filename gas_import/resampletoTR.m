%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [resamp_time,resamp_data] = resampletoTR(TR,time,data)
% this function resamples the "data" based on the "TR" given as input and
% spits out the resampled time points and the corresponding data

ntime = time - time(1,1); %make the first time value correspond to 0
start_time = ntime(1,1);
end_time = ntime(end,1);

%resample the data
resamp_time = linspace(start_time,end_time,floor(end_time/TR));
resamp_data = interp1(ntime,data,resamp_time); 
end