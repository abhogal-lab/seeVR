function quickPlot(vol, nRows, nCols, cmap)
% quickPlot  —  fast visual check of a 3-D volume
%
%   quickPlot(vol)                           → 4×4 grid, auto-colormap
%   quickPlot(vol, nRows, nCols)             → custom grid
%   quickPlot(vol, nRows, nCols, cmap)       → explicit colormap
%
%   • Only the 3rd dimension is treated as “slice”.
%   • Blank slices and blank XY borders are skipped/cropped.
%   • Each subplot autoscale uses imagesc(volSlice, []).

% ----------------------------------------------------------------------
% Alex Bhogal  ·  20-Jun-2025
% ----------------------------------------------------------------------

% -------------- defaults & sanity checks -------------------------------
if nargin < 2 || isempty(nRows), nRows = 4; end
if nargin < 3 || isempty(nCols), nCols = 4; end
assert(ndims(vol) == 3, 'Input must be a 3-D matrix.');
vol(vol == 0) = NaN; %for black background  regardless of colormap. If you need zeros to be show, comment out
% -------------- automatic colormap choice ------------------------------
if nargin < 4 || isempty(cmap)
    if islogical(vol) || numel(unique(vol(~isnan(vol)))) <= 3
        cmap = gray;
    else
        cmap = parula;
    end
elseif ischar(cmap) || isstring(cmap)
    cmap = feval(char(cmap),256);
end

vol = double(vol);

% -------------- find informative slice range ---------------------------
sliceSum = squeeze(nansum(nansum(abs(vol), 1), 2));
nzSlices = find(sliceSum > 0);

if isempty(nzSlices)
    warning('quickPlot:NoData','Volume contains no non-zero voxels.');
    return
end

firstSlice = nzSlices(1);
lastSlice  = nzSlices(end);

% -------------- choose slices to display -------------------------------
maxPlots = nRows*nCols;
if numel(nzSlices) <= maxPlots
    sliceIdx = nzSlices;
else
    sliceIdx = round(linspace(firstSlice,lastSlice,maxPlots));
end
maxPlots = numel(sliceIdx);

% -------------- optional XY bounding-box crop --------------------------
maskXY = any(vol(:,:,sliceIdx) ~= 0 & ~isnan(vol(:,:,sliceIdx)), 3);
[rows, cols] = find(maskXY);
if isempty(rows)
    r1=1; r2=size(vol,1); c1=1; c2=size(vol,2);
else
    pad = 4;
    r1 = max(min(rows)-pad,1);  r2 = min(max(rows)+pad,size(vol,1));
    c1 = max(min(cols)-pad,1);  c2 = min(max(cols)+pad,size(vol,2));
end

% -------------- plot ---------------------------------------------------
fig = figure('Color','k', ...
             'Name','quickPlot', ...
             'Units','pixels', ...
             'Position',[100 100 1400 800]);  % ← Wider figure

% main tiled layout for slices
t = tiledlayout(nRows, nCols, 'Padding','compact', 'TileSpacing','compact');

% plot each selected slice
for k = 1:maxPlots
    ax = nexttile;
    slice = vol(r1:r2, c1:c2, sliceIdx(k));
    img = imagesc(ax, slice);
    set(img, 'AlphaData', ~isnan(slice));  % hide NaNs
    axis(ax,'image','off');
    colormap(ax, cmap);
    set(ax, 'Color', 'k');                % black background
    title(ax, sprintf('z = %d',sliceIdx(k)), 'Color','w', 'FontSize',8);
end

% attach colorbar on side using dummy axis
cbAx = axes('Parent',fig,'Position',[0.93 0.15 0.02 0.7],'Visible','off');
colormap(cbAx,cmap);
caxis(cbAx,[min(vol(:),[],'omitnan') max(vol(:),[],'omitnan')]);
cb = colorbar(cbAx);
cb.Units = 'normalized';
cb.Position = [0.94 0.15 0.015 0.7];
cb.Color = 'w';
cb.FontSize = 10;
cb.Label.String = 'Intensity';

% sgtitle('quickPlot – fast slice overview', 'Color','w');
end
