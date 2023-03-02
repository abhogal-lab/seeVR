% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <resampleGEN3: resamples GEN3 breathing data to the time-step defined by TR >
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

function[BBBdata,yi yyi] = resampleGEN3(seqpath,TR,dyn,oversample)

cd(seqpath); 
traceFile =  ls('*bbb_mod*')

%import BBB
fileToRead = [seqpath,traceFile];
DELIMITER = '\t';
HEADERLINES = 1;
BBBdata = import_respfiles(fileToRead, DELIMITER, HEADERLINES);
petCO2 = BBBdata.data(:,2); petO2 = BBBdata.data(:,3); time = BBBdata.data(:,1);
time = time-time(1,1);

xdata = [TR:TR:TR*dyn];
resample_rate = TR*oversample; 
start_time = 60*time(1,1);
end_time = 60*time(end,1);
total_time = end_time - start_time;
    y = petCO2; %respiract raw CO2 trace
    yy = petO2; %respiract raw CO2 trace
    x = 60*time; %respiract time vector
   % xi = linspace(start_time,(floor(total_time/resample_rate)-1)*resample_rate,floor(total_time/resample_rate)); %respiract time vector
     xi = linspace(start_time,(floor(total_time/resample_rate))*resample_rate,floor(total_time/resample_rate)); %respiract time vector
    yi = interp1(x,y,xi); %yi is resampled CO2 
    yyi = interp1(x,yy,xi); %yi is resampled O2

end