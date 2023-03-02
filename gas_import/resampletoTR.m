% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <resampletoTR: resamples breathing data to TR >
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

function [resamp_time,resamp_data] = resampletoTR(TR,time,data)
% this function resamples the "data" based on the "TR" given as input and
% spits out the resampled time points and the corresponding data
if iscolumn(time); else; time = time'; end
if iscolumn(data); else; data = data'; end

ntime = time - time(1,1); %make the first time value correspond to 0
start_time = ntime(1,1);
end_time = ntime(end,1);

%resample the data
resamp_time = linspace(0,end_time,floor(end_time/TR)+1);
resamp_data = interp1(ntime,data,resamp_time); 
end