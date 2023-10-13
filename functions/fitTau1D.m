% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <fitTau1D: calculates the dispersion time constant between in input
% signal and an associated MRI signal timeseries
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
%
% probe: usually PetCO2/O2 trace
% 
% ts: input timeseries data (i.e. average BOLD signal calculated over an ROI using meanTimeseries)
%

function [y_fit, b] = fitTau1D(probe, ts, opts)
global opts
ts = double(ts);
if iscolumn(probe); probe = probe'; end
if iscolumn(ts); ts = ts'; end

if isfield(opts,'maxTau'); else; opts.maxTau = 300; end                   %maximum exponential dispersion time constant - data dependent

opts.dynamicdir = fullfile(opts.resultsdir,'tau'); mkdir(opts.dynamicdir);

t = double(opts.TR:opts.TR:length(probe)*opts.TR);


options = optimoptions('lsqcurvefit','Display','none','FunctionTolerance',1.0000e-8,...
    'StepTolerance', 1.0000e-8, 'MaxIter',150);

nr_params = 3;

b = nan([length(probe), nr_params]);

model = (@(a,t) a(1)*rescale(real(ifft(ifftshift(fftinput(probe).*fftexponential(a(2),t)))))+a(nr_params));
a0 = [0 20 0];
lb = [ 0  0 -Inf];
ub = [ Inf opts.maxTau Inf];

tic
try
    b = lsqcurvefit(model,[range(ts) 20 0],t,ts,lb,ub,options)
catch
    %continue
    disp('error; please check inputs')
end
y_fit = model(b,t);
% Plot the results
figure
plot(t,(ts),'r-',t,(y_fit),'b-');
legend('Observed Data','Fitted Model');

end


