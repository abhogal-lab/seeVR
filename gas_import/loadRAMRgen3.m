%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [corrvec_CO2 corrvec_O2] = loadRAMRgen3(opts)
%written by Alex Bhogal a.bhogal@umcutrecht.nl
%functions to load respiract endtidal traces from 3rd generation respiract
%system. Input is path to directory where gas traces are saved. Output is
%resampled O2 and CO2 traces at the desired sampling frequency (TR of scan)

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

