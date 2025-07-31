function quickPlot(vol, nRows, nCols, scale, cmap, nzSlices)
% quickPlot  —  fast visual check of a 3‑D volume                  (v1.4)
%
%   quickPlot(vol)                            → 4×4 grid, auto‑colormap
%   quickPlot(vol, nRows, nCols)              → custom grid
%   quickPlot(vol, nRows, nCols, scale)       → explicit colour scale
%   quickPlot(vol, nRows, nCols, scale, cmap) → explicit scale + map
%   quickPlot(vol, ..., nzSlices)             → manually define slices to show
%
%   •  The 3rd dimension is treated as “slice”.
%   •  Blank slices and blank XY borders are skipped/cropped.
%   •  If  scale = []  (or is omitted) each slice autoscales exactly
%      like  imagesc(slice,[])  and the colour‑bar spans
%         [min(vol(:))  max(vol(:))].
%   •  If  scale  is provided and **symmetric about zero** (scalar
%      shorthand allowed), every subplot and the shared colour‑bar use
%      those limits so the physical mid‑point corresponds to 0, and a
%      tick label “0” is added at the centre.

% ---------- defaults & sanity checks -----------------------------------
if nargin < 2 || isempty(nRows),  nRows = 4;   end
if nargin < 3 || isempty(nCols),  nCols = 4;   end
if nargin < 4,                    scale = [];  end
if nargin < 5,                    cmap  = [];  end
if nargin < 6,                    nzSlices = []; end

assert(ndims(vol) == 3, 'quickPlot:Input', 'Input must be a 3‑D matrix.');

vol(vol == 0) = NaN;
vol          = double(vol);

% ---------- automatic / explicit colormap ------------------------------
if isempty(cmap)
    if islogical(vol) || numel(unique(vol(~isnan(vol)))) <= 3
        cmap = gray(256);
    else
        cmap = parula(256);
    end
elseif ischar(cmap) || isstring(cmap)
    cmap = feval(char(cmap), 256);
end

% ---------- decide what to do with SCALE -------------------------------
symScale  = false;
userScale = ~isempty(scale);

if userScale
    if numel(scale) == 1
        scale     = [-abs(scale) abs(scale)];
        symScale  = true;
    elseif numel(scale) == 2 && abs(scale(1) + scale(2)) < eps
        scale     = scale(:)';
        symScale  = true;
    else
        scale     = scale(:)';
    end
end

% ---------- find or use user-specified slice range ---------------------
if isempty(nzSlices)
    sliceSum = squeeze(nansum(nansum(abs(vol), 1), 2));
    nzSlices = find(sliceSum > 0);

    if isempty(nzSlices)
        warning('quickPlot:NoData', 'Volume contains no non‑zero voxels.');
        return
    end

    fprintf('Non-zero slice range: z = %d to %d\n', nzSlices(1), nzSlices(end));
else
    if isvector(nzSlices) && numel(nzSlices) == 2
        nzSlices = nzSlices(1):nzSlices(2);
    end
    nzSlices = unique(nzSlices(nzSlices >= 1 & nzSlices <= size(vol,3)));

    if isempty(nzSlices)
        warning('quickPlot:InvalidInput', 'Provided nzSlices are out of volume bounds.');
        return
    end

    fprintf('Using user-specified slice range: z = %d to %d\n', nzSlices(1), nzSlices(end));
end

firstSlice = nzSlices(1);
lastSlice  = nzSlices(end);

% ---------- choose slices to display -----------------------------------
maxPlots = nRows * nCols;
if numel(nzSlices) <= maxPlots
    sliceIdx = nzSlices;
else
    sliceIdx = round(linspace(firstSlice, lastSlice, maxPlots));
end
maxPlots = numel(sliceIdx);

% ---------- optional XY bounding box crop ------------------------------
maskXY       = any(vol(:, :, sliceIdx) ~= 0 & ~isnan(vol(:, :, sliceIdx)), 3);
[rows, cols] = find(maskXY);
if isempty(rows)
    r1 = 1;               r2 = size(vol, 1);
    c1 = 1;               c2 = size(vol, 2);
else
    pad = 4;
    r1 = max(min(rows) - pad, 1);               r2 = min(max(rows) + pad, size(vol, 1));
    c1 = max(min(cols) - pad, 1);               c2 = min(max(cols) + pad, size(vol, 2));
end

% ----------------------------- plot ------------------------------------
fig = figure('Color', 'k', 'Name', 'quickPlot', 'Units', 'pixels', ...
             'Position', [100 100 1400 800]);

t = tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:maxPlots
    ax    = nexttile;
    slice = vol(r1:r2, c1:c2, sliceIdx(k));

    img   = imagesc(ax, slice);
    set(img, 'AlphaData', ~isnan(slice));
    axis(ax, 'image', 'off');
    colormap(ax, cmap);

    if userScale
        caxis(ax, scale);
    end

    title(ax, sprintf('z = %d', sliceIdx(k)), ...
              'Color', 'w', 'FontSize', 8);
end

% -------------------------- shared colour‑bar --------------------------
cbAx = axes('Parent', fig, 'Position', [0.93 0.15 0.02 0.7], 'Visible', 'off');
colormap(cbAx, cmap);

if userScale
    caxis(cbAx, scale);
else
    caxis(cbAx, [min(vol(:), [], 'omitnan')  max(vol(:), [], 'omitnan')]);
end

cb = colorbar(cbAx);
cb.Units        = 'normalized';
cb.Position     = [0.94 0.15 0.015 0.7];
cb.Color        = 'w';
cb.FontSize     = 10;
cb.Label.String = 'Intensity';

if symScale
    cb.Ticks = [scale(1) 0 scale(2)];
end

% sgtitle('quickPlot – fast slice overview', 'Color', 'w');  % optional
end
