function quickPlot(vol, nRows, nCols)
% quickPlot  —  fast visual check of a 3-D volume
%
%   quickPlot(vol)              shows a 4×4 grid of slices
%   quickPlot(vol, nRows, nCols) lets you choose the grid size
%
%   • Only the 3rd dimension is interpreted as “slice”.
%   • Slices that are entirely zero or NaN are skipped automatically.
%   • Each subplot is autoscaled (imagesc with []) and shown in gray.
%
% ----------------------------------------------------------------------
% Alex Bhogal  ·  20-Jun-2025
% ----------------------------------------------------------------------

% ---------------- defaults & sanity checks ----------------------------
if nargin < 2 || isempty(nRows), nRows = 4; end
if nargin < 3 || isempty(nCols), nCols = 4; end
assert(ndims(vol) == 3, 'Input must be a 3-D matrix.');

vol = double(vol);                       % ensure numeric

% ---------------- find informative slice range ------------------------
sliceSum = squeeze(sum(sum(abs(vol),1),2));     % energy per slice
nzSlices = find(sliceSum > 0 & ~isnan(sliceSum));

if isempty(nzSlices)
    warning('quickPlot:NoData','Volume contains no non-zero voxels.');
    return
end

firstSlice = nzSlices(1);
lastSlice  = nzSlices(end);

% ---------------- choose slices to display ----------------------------
nPlots  = nRows * nCols;
if numel(nzSlices) <= nPlots
    sliceIdx = nzSlices;                       % show all non-zero slices
else
    sliceIdx = round(linspace(firstSlice, lastSlice, nPlots));
end
nPlots = numel(sliceIdx);                      % may be fewer than 16

% ---------------- optional XY bounding-box crop -----------------------
% (crops empty borders, keeps at least 4-voxel padding)
maskXY = any(vol(:,:,sliceIdx) ~= 0 & ~isnan(vol(:,:,sliceIdx)), 3);
[rows, cols] = find(maskXY);
if isempty(rows)
    r1 = 1; r2 = size(vol,1);
    c1 = 1; c2 = size(vol,2);
else
    pad  = 4;
    r1 = max(min(rows)-pad,1);  r2 = min(max(rows)+pad, size(vol,1));
    c1 = max(min(cols)-pad,1);  c2 = min(max(cols)+pad, size(vol,2));
end

% ------------------- plot ---------------------------------------------
figure('Color','k','Name','quickPlot');
for k = 1:nPlots
    subplot(nRows, nCols, k);
    imagesc(vol(r1:r2, c1:c2, sliceIdx(k)));   %#ok<NASGU>
    axis image off;
    colormap gray;
    title(sprintf('z = %d', sliceIdx(k)),'Color',[1 1 1],'FontSize',8);
end
sgtitle('quickPlot – fast slice overview','Color',[1 1 1]);

end
