% Copyright (C) Alex A. Bhogal, 2025
% a.bhogal@umcutrecht.nl
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
function quickPlot_overlay(volBase, volOverlay, nRows, nCols, ...
                           scaleBase, cmapBase, scaleOverlay, cmapOverlay, ...
                           alphaOverlay, nzSlices)
% quickPlot_overlay — overlay visualization for 3D volumes (no panel box)
%
% Two colorbars (Base, Overlay) now occupy 80% height, clean black bg.

% ---------- defaults ----------------------------------------------------
if nargin < 3 || isempty(nRows),         nRows = 4;            end
if nargin < 4 || isempty(nCols),         nCols = 4;            end
if nargin < 5,                           scaleBase = [];       end
if nargin < 6,                           cmapBase  = [];       end
if nargin < 7,                           scaleOverlay = [];    end
if nargin < 8,                           cmapOverlay  = [];    end
if nargin < 9 || isempty(alphaOverlay),  alphaOverlay = 0.5;   end
if nargin < 10,                          nzSlices = [];        end

assert(ndims(volBase)==3, 'quickPlot_overlay:Input', 'Base input must be a 3-D matrix.');
if ~isempty(volOverlay)
    assert(ndims(volOverlay)==3, 'quickPlot_overlay:Input', 'Overlay input must be a 3-D matrix.');
    assert(isequal(size(volBase), size(volOverlay)), ...
        'quickPlot_overlay:SizeMismatch','Base and overlay volumes must have the same size.');
end

volBase = double(volBase);     volBase(volBase==0) = NaN;
hasOverlay = ~isempty(volOverlay) && any(~isnan(volOverlay(:)));
if hasOverlay
    volOverlay = double(volOverlay); volOverlay(volOverlay==0) = NaN;
end

% ---------- colormaps & scales ------------------------------------------
cmapBase    = localResolveColormap(cmapBase, volBase);
if hasOverlay
    cmapOverlay = localResolveColormap(cmapOverlay, volOverlay);
end
[userScaleB, symScaleB, scaleBase]       = localResolveScale(scaleBase);
[userScaleO, symScaleO, scaleOverlay]     = localResolveScale(scaleOverlay);

% ---------- slice range -------------------------------------------------
if isempty(nzSlices)
    sliceSum = squeeze(sum(sum(abs(volBase),1),2));
    if hasOverlay, sliceSum = sliceSum + squeeze(sum(sum(abs(volOverlay),1),2)); end
    nzSlices = find(sliceSum > 0);
else
    if isvector(nzSlices) && numel(nzSlices)==2, nzSlices = nzSlices(1):nzSlices(2); end
end
if isempty(nzSlices)
    warning('quickPlot_overlay:NoData','Volumes contain no non-zero voxels.'); return;
end
sliceIdx = round(linspace(nzSlices(1), nzSlices(end), min(nRows*nCols,numel(nzSlices))));

% ---------- crop region -------------------------------------------------
if hasOverlay
    maskXY = any((~isnan(volBase(:,:,sliceIdx)) & volBase(:,:,sliceIdx)~=0) | ...
                 (~isnan(volOverlay(:,:,sliceIdx)) & volOverlay(:,:,sliceIdx)~=0), 3);
else
    maskXY = any(~isnan(volBase(:,:,sliceIdx)) & volBase(:,:,sliceIdx)~=0, 3);
end
[rows, cols] = find(maskXY);
if isempty(rows), r1=1; r2=size(volBase,1); c1=1; c2=size(volBase,2);
else, pad=4;
    r1=max(min(rows)-pad,1); r2=min(max(rows)+pad,size(volBase,1));
    c1=max(min(cols)-pad,1); c2=min(max(cols)+pad,size(volBase,2));
end

% ---------- figure & layout (extra right margin for colorbars + labels) ----------
fig = figure('Color','k', 'InvertHardCopy','off', ...   % <— added
             'Name','quickPlot_overlay','Units','pixels', ...
             'Position',[100 100 800 800]);

% Tiled layout in a narrower panel to leave ample right margin
tiledPanel = uipanel('Parent', fig, 'Units', 'normalized', ...
                     'Position', [0.04 0.05 0.76 0.90], ...   % 80% width = more right space
                     'BackgroundColor','k', 'BorderType','none');

t = tiledlayout(tiledPanel, nRows, nCols, ...
                'Padding','compact', 'TileSpacing','compact');

% ---------- render tiles ------------------------------------------------
for k = 1:numel(sliceIdx)
    ax = nexttile(t);
    set(ax,'Color','k');    
    baseSlice = volBase(r1:r2, c1:c2, sliceIdx(k));
    [RGBb,~,~] = localMapToRGB(baseSlice, cmapBase, scaleBase, userScaleB);

    if hasOverlay
        ovSlice = volOverlay(r1:r2, c1:c2, sliceIdx(k));
        [RGBo,~,~] = localMapToRGB(ovSlice, cmapOverlay, scaleOverlay, userScaleO);
        A = alphaOverlay * double(~isnan(ovSlice));
        out = RGBb;
        for ch = 1:3
            out(:,:,ch) = (1 - A).*RGBb(:,:,ch) + A.*RGBo(:,:,ch);
        end
        image(ax, out);
    else
        image(ax, RGBb);
    end
    axis(ax,'image','off');
    title(ax, sprintf('z = %d', sliceIdx(k)), 'Color','w','FontSize',8);

end

% ---------- dual colour strips with labels (no ticks/box) ----------
barWidth   = 0.014;
barHeight  = 0.80;
barBottom  = 0.10;

marginFromTiles = 0.025;   % bars further LEFT when larger
gap             = 0.060;   % spacing between the two bars

panelRight = tiledPanel.Position(1) + tiledPanel.Position(3);
xBase      = panelRight + marginFromTiles;         % base strip (left)
xOverlay   = xBase + barWidth + gap;               % overlay strip (right)

% safe limits (like localSafeCaxis)
getLims = @(s,v) ( ...
    ~isempty(s) && numel(s)==2 && all(isfinite(s)) ) .* double(s) + ...
    ( isempty(s) || numel(s)~=2 || ~all(isfinite(s)) ) .* ...
    [min(v(:),[],'omitnan') max(v(:),[],'omitnan')];

limsB = getLims(scaleBase,    volBase);  if limsB(2) <= limsB(1), limsB = limsB + [-1 1]*1e-6; end
if hasOverlay
    limsO = getLims(scaleOverlay, volOverlay);     if limsO(2) <= limsO(1), limsO = limsO + [-1 1]*1e-6; end
end

% draw function: gradient strip + right-side labels (no ticks, no box)
function drawStrip(xpos, cmap, lims, symFlag)
    ax = axes('Parent',fig, 'Position',[xpos barBottom barWidth barHeight], ...
              'Color','k','Box','off');                   % no white box
    vals = linspace(lims(1), lims(2), 256).';
    imagesc(ax, [0 1], [lims(1) lims(2)], [vals vals]);   % vertical strip
    set(ax,'YDir','normal');                              % low→bottom
    colormap(ax, cmap);  caxis(ax, lims);
    axis(ax, 'off');                                      % hide ticks & axes

    % pick label positions/format
    if symFlag, labelVals = [lims(1) 0 lims(2)];
    else,       labelVals = linspace(lims(1), lims(2), 6);
    end
    step = mean(diff(labelVals));
    if ~isfinite(step), fmt = '%.3g';
    elseif abs(step) < 1,  fmt = '%.2f';
    elseif abs(step) < 10, fmt = '%.1f';
    else,                  fmt = '%.0f';
    end

    labelOffset = 0.08;  % move labels further right (normalized units)
    hold(ax,'on');
    for v = labelVals
        y = (v - lims(1)) / max(eps, (lims(2)-lims(1)));   % normalize y
        text(ax, 1 + labelOffset, y, sprintf(fmt, v), ...
             'Units','normalized', 'Color','w', ...
             'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
             'Clipping','off', 'FontSize',10);
    end
    hold(ax,'off');
end

% Base strip + labels
drawStrip(xBase,    cmapBase,    limsB, symScaleB);

% Overlay strip + labels
if hasOverlay
    [~, symO, ~] = localResolveScale(scaleOverlay);
    drawStrip(xOverlay, cmapOverlay, limsO, symO);
end



end
% ========================== helpers ====================================
function cmap = localResolveColormap(cmap, vol)
    if isempty(cmap)
        u = unique(vol(~isnan(vol)));
        if islogical(vol) || numel(u)<=3, cmap=gray(256);
        else, cmap=parula(256); end
    elseif ischar(cmap) || isstring(cmap)
        cmap = feval(char(cmap),256);
    elseif isnumeric(cmap)
        if size(cmap,2)~=3, error('Colormap must be N×3'); end
        if size(cmap,1)~=256
            xi=1:256; x=linspace(1,256,size(cmap,1));
            cmap=[interp1(x,cmap(:,1),xi)',interp1(x,cmap(:,2),xi)',interp1(x,cmap(:,3),xi)'];
            cmap=max(min(cmap,1),0);
        end
    end
end

function [userScale,symScale,scaleOut]=localResolveScale(scaleIn)
    userScale=~isempty(scaleIn); symScale=false; scaleOut=scaleIn;
    if ~userScale, return; end
    if numel(scaleIn)==1
        scaleOut=[-abs(scaleIn) abs(scaleIn)]; symScale=true;
    elseif numel(scaleIn)==2
        scaleOut=double(scaleIn(:))';
        if scaleOut(2)<scaleOut(1), scaleOut=fliplr(scaleOut); end
        if scaleOut(1)==scaleOut(2), scaleOut=scaleOut+[-1 1]*1e-6; end
        if abs(scaleOut(1)+scaleOut(2))<eps, symScale=true; end
    else
        error('Scale must be [] (auto), scalar, or [min max]');
    end
end

function [RGB,cmin,cmax]=localMapToRGB(slice,cmap,scale,userScale)
    if all(isnan(slice(:)))
        RGB = zeros([size(slice) 3],'like',double(slice));
        cmin = 0; cmax = 1; 
        return
    end

    nanMask = isnan(slice);
    if userScale && ~isempty(scale)
        cmin=scale(1); cmax=scale(2);
    else
        cmin=min(slice(:),[],'omitnan'); 
        cmax=max(slice(:),[],'omitnan');
        if cmin==cmax, cmin=cmin-1; cmax=cmax+1; end
    end

    f = (slice - cmin) / (cmax - cmin);
    f = max(min(f,1),0);

    idx = 1 + floor(f*255);
    idx(nanMask) = 1;                % temp index (we’ll overwrite with black)
    RGB = ind2rgb(idx, cmap);

    % Force NaN pixels to black (don’t rely on cmap(1,:))
    if any(nanMask(:))
        m3 = repmat(nanMask, 1, 1, 3);
        RGB(m3) = 0;
    end
end


function localSafeCaxis(ax,scale,vol)
    if ~isempty(scale)&&numel(scale)==2&&isfinite(scale(1))&&isfinite(scale(2))
        lo=scale(1);hi=scale(2); if hi<=lo, [lo,hi]=deal(lo-1,hi+1); end
    else
        lo=min(vol(:),[],'omitnan'); hi=max(vol(:),[],'omitnan');
        if ~isfinite(lo)||~isfinite(hi)||hi<=lo, lo=-1; hi=1; end
    end
    caxis(ax,[lo hi]);
end
