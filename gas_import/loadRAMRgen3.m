% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <loadRAMRgen3: imports GEN3 RespirAct data >
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
% This functions loads respiract endtidal traces from 3rd generation respiract
% system. 
% 
% opts.seqpath: sequence path pointing to saved GEN3 respiratory data
%
% opts.TR: TR of MRI data - traces are resampled to match
%
% opts.dyn: number of datapoints in MRI timeseries data
%
% corrvec_CO2: the temporally resampled (to opts.TR) end-expired CO2 traces
%
% corrvec_O2: the temporally resampled (to opts.TR) end-expired O2 traces
function [corrvec_CO2 corrvec_O2] = loadRAMRgen3(opts)


cd(opts.seqpath)
oversample = 1;
[BBBdata corrvec_CO2 corrvec_O2] = resampleGEN3(opts.seqpath,opts.TR,opts.dyn,oversample)
%import RAW
rawFile = ls('*raw*')
raw = importRAW(rawFile);

%plot time versus BBB, raw is: 1)time, 2)Pmouth, 3)Pdelay, 4)PCO2, 5)PO2
figure;
sz = 140;
plot(table2array(raw(:,1)),table2array(raw(:,4)),'k', 'LineWidth', 0.5); hold on;
scatter(BBBdata.data(:,1),BBBdata.data(:,2),sz, 'm.');
figure; %O2
sz = 140;
plot(table2array(raw(:,1)),table2array(raw(:,5)),'k', 'LineWidth', 0.5); hold on
scatter(BBBdata.data(:,1),BBBdata.data(:,3),sz, 'c.');
saveas(gcf,[opts.figdir,'endtidal_data.fig']);
end

