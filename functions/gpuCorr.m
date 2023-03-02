function [data] = gpuCorr(a,b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
a2 = gpuArray(a);
c2 = gpuArray(repmat(b,1, size(a,1)));
data = xcorr2(a2,c2);
end

