% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <prepNuisance: adds nuisance signal derivatives, drift term, and input probe dispersion>
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
% This function perpares nuisance signals before the data scrubbing process.
%
% nuisance: an array of nuisance regressors (i.e. motion,
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
function [keep, leave] = prepNuisance(nuisance,probe, opts)
global opts

if isfield(opts,'filter_regressors'); else; opts.filter_regressors = 1; end
if isfield(opts,'add_drift'); else; opts.add_drift = 1; end % add a linear drift term
if isfield(opts,'legendre_order'); else; opts.legendre_order = 1; end % increasing this adds additional legendre polynomials
if isfield(opts,'add_derivatives'); else; opts.add_derivatives = 2; end % if = 1 motion trans/rot, if = 2 or 3 then add time derivative and square, respectively
if isfield(opts,'motioncorr'); else; opts.motioncorr  = 0.3; end %exclude traces with higher correlation values
if opts.motioncorr > 1 || opts.motioncorr < 0
    opts.motioncorr = 0.3;
end

% generate derivatives
if opts.add_derivatives > 4 || opts.add_derivatives < 1
    opts.add_derivatives = 1;
end

switch opts.add_derivatives
    case 1
        disp('using input nuisance signals')
        motion = nuisance;
    case 2
        disp('using input nuisance signals and temporal derivative')
        dtnuisance =  gradient(nuisance); % temporal derivative
        motion = [nuisance dtnuisance];
    case 3
        disp('using input nuisance signals and temporal derivative and square of input signal')
        sqnuisance = nuisance.*nuisance; % square of motion
        dtnuisance =  gradient(nuisance); % temporal derivative
        motion = [nuisance dtnuisance sqnuisance];
end

if opts.add_drift ~= 0
    disp('adding drift term')
    if opts.legendre_order == 0; opts.legendre_order = 1; end
    for ii=1:opts.legendre_order
    drift_term(ii,:) = LegendreN(ii, 1:1:length(nuisance));
    end
    %concatenate motion params with drift term
    %np0 = [motion drift_term' flip(drift_term')];
    np0 = [motion drift_term'];

else
    disp('no drift term added')
    np0 = motion;
end

%rescale nuisance parameters
disp('rescaling nuisance between -1 & 1')
for kk=1:size(np0,2)
    np0(:,kk) = rescale(np0(:,kk),-1,1);
end

% filter regressors
if opts.filter_regressors
    disp('filtering regressors based on correlation value supplied in opts.motioncorr')
    [keep, leave] = filtRegressor(np0, probe, opts);
end
end

