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