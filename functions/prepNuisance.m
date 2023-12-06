function [keep] = prepNuisance(nuisance,reference, opts)
global opts

if isfield(opts,'filloutliers'); else; opts.filloutliers = 1; end
if isfield(opts,'filter_regressors'); else; opts.filter_regressors = 1; end
if isfield(opts,'add_drift'); else; opts.add_drift = 1; end % add a linear drift term
if isfield(opts,'legendre_order'); else; opts.legendre_order = 1; end % increasing this adds additional legendre polynomials
if isfield(opts,'add_derivatives'); else; opts.add_derivatives = 2; end % if = 1 motion trans/rot, if = 2 or 3 then add time derivative and square, respectively
if isfield(opts,'motioncorr'); else; opts.motioncorr  = 0.4; end %exclude traces with higher correlation values
if opts.motioncorr > 1 || opts.motioncorr < 0
    opts.motioncorr = 0.4;
end

%check for nonesense case
if opts.add_derivatives == 0 && opts.add_drift == 0
    disp('options indicate no nuisance signals should be included')
    disp('this is probably in error - check options')
    disp('for this run, rotations and translations will be included')
    disp('opts.add_derivatives = 1')
    opts.add_derivatives = 1;
end

% remove outliers
if opts.filloutliers == 1
    try
        for kk=1:size(nuisance,2)
            nuisance(:,kk) = filloutliers(nuisance(:,kk),'center');
        end
        
    catch
        disp('error, skipping outlier removal')
    end
end

% generate derivatives
np = [];
dtnuisance =  gradient(nuisance); % temporal derivative
sqnuisance = nuisance.*nuisance; % square of motion

if opts.add_derivatives > 4 || opts.add_derivatives < 1
    opts.add_derivatives = 0;
end

switch opts.add_derivatives
    case 0
        motion = [];
    case 1
        motion = nuisance;
    case 2
        motion = [nuisance dtnuisance];
    case 3
        motion = [nuisance dtnuisance sqnuisance];
end

if opts.add_drift
    drift_term = LegendreN(opts.legendre_order, 1:1:length(nuisance));
    %concatenate motion params with drift term
    np0 = [motion drift_term'];
else
    np0 = motion;
end

%demean and rescale nuisance parameters
for kk=1:size(np0,2)
    np0(:,kk) = rescale(demeanData(np0(:,kk)),-1,1);
end

% filter regressors
if opts.filter_regressors
[keep, ~] = filtRegressor(np0, reference, opts);
end
end

