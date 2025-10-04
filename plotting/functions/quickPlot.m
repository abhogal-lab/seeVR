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
    sliceSum = squeeze(sum(sum(abs(v),1),2));
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
fig = figure('Color','k','InvertHardCopy','off', ...
             'Name','quickPlot','Units','pixels', ...
             'Position',[100 100 800 800]);

t = tiledlayout(fig, nRows, nCols, 'Padding','compact', 'TileSpacing','compact');

% Reserve a right gutter for the colorbar (left x, bottom y, width, height)
% ↓ Reduce the width to move the image grid left (e.g., 0.76–0.80 works well)
t.OuterPosition = [0.04 0.05 0.78 0.90];

for k = 1:maxPlots
    ax    = nexttile(t);
    slice = vol(r1:r2, c1:c2, sliceIdx(k));
    img   = imagesc(ax, slice);
    set(img, 'AlphaData', ~isnan(slice));
    axis(ax, 'image', 'off');
    colormap(ax, cmap);
    if userScale, caxis(ax, scale); end
    title(ax, sprintf('z = %d', sliceIdx(k)), 'Color','w','FontSize',8);
end

% -------------------------- shared colour-bar (labels only) -------------
gutter     = 0.015;      % gap from tile area to bar (smaller → more left)
barBottom  = 0.10;
barHeight  = 0.80;
barWidth   = 0.018;

panelRight = t.OuterPosition(1) + t.OuterPosition(3);   % right edge of tiles
xBar       = panelRight + gutter;

% Colour limits
if userScale
    clim = scale;
else
    clim = [min(vol(:),[],'omitnan')  max(vol(:),[],'omitnan')];
end
if ~isfinite(clim(1)) || ~isfinite(clim(2)) || clim(2) <= clim(1)
    clim = [0 1];
end

% --- gradient strip as the "bar" ---
cbAX = axes('Parent',fig, 'Position',[xBar barBottom barWidth barHeight], ...
            'Color','k','Box','off');
vals = linspace(clim(1),clim(2),256).';
imagesc(cbAX, [0 1], [clim(1) clim(2)], [vals vals]);   % vertical strip
set(cbAX,'YDir','normal');           % low at bottom, high at top
colormap(cbAX, cmap);
caxis(cbAX, clim);

% Hide ticks/axes completely (no tick marks)
axis(cbAX,'off');

% --- draw labels to the RIGHT (no ticks) ---
if symScale
    labelVals = [clim(1) 0 clim(2)];
else
    labelVals = linspace(clim(1),clim(2),6);  % 6 labels by default
end

% choose a sensible numeric format
step = mean(diff(labelVals));
if step < 1,   fmt = '%.2f';
elseif step < 10, fmt = '%.1f';
else            fmt = '%.0f';
end

labelOffset = 0.2;   % ↑ move labels further right (normalized X units)
hold(cbAX,'on');
for v = labelVals
    y = (v - clim(1)) / max(eps, (clim(2) - clim(1)));  % normalize to [0,1]
    text(cbAX, 1 + labelOffset, y, sprintf(fmt, v), ...
         'Units','normalized', 'Color','w', ...
         'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
         'Clipping','off', 'FontSize',10);
end
hold(cbAX,'off');


%cb.Label.Position = [3 0 0];   % nudge label further to the right

% sgtitle('quickPlot – fast slice overview', 'Color', 'w');  % optional
end