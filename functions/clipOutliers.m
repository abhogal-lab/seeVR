function vol = clipOutliers(vol, prcLo, prcHi)
% Replace values below/above the chosen percentiles with the
% percentile itself.  Keeps volume size and NaNs intact.
mask          = isfinite(vol);           % ignore NaN / Inf
vec           = vol(mask);
limits        = prctile(vec, [prcLo prcHi]);   % e.g. [0.5 99.5]
vol(mask & vol<limits(1)) = limits(1);
vol(mask & vol>limits(2)) = limits(2);
end