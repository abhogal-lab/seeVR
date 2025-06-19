% Copyright (C) Alex A. Bhogal, 2025, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <quickplot: quick and dirty way to look at 3dvol data>
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
%
% ----------------------------------------------------------------------
% Alex Bhogal  ·  20-Jun-2025
% ----------------------------------------------------------------------

% -------------- defaults & sanity checks -------------------------------
if nargin < 2 || isempty(nRows), nRows = 4; end
if nargin < 3 || isempty(nCols), nCols = 4; end
assert(ndims(vol) == 3, 'Input must be a 3-D matrix.');

% -------------- automatic colormap choice ------------------------------
if nargin < 4 || isempty(cmap)
    % Heuristic “is this a mask?”
    if islogical(vol) || numel(unique(vol(~isnan(vol)))) <= 3
        cmap = gray;                         % mask → gray
    else
        cmap = parula;                       % normal data → parula
    end
elseif ischar(cmap) || isstring(cmap)
    cmap = feval(char(cmap),256);            % convert name to array
end

vol = double(vol);

% -------------- find informative slice range ---------------------------
sliceSum = squeeze(sum(sum(abs(vol),1),2));
nzSlices = find(sliceSum > 0 & ~isnan(sliceSum));

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
figure('Color','k','Name','quickPlot');
for k = 1:maxPlots
    subplot(nRows,nCols,k);
    imagesc(vol(r1:r2,c1:c2,sliceIdx(k))); 
    axis image off; 
    colormap(cmap); 
    title(sprintf('z = %d',sliceIdx(k)),'Color',[1 1 1],'FontSize',8);
end
sgtitle('quickPlot – fast slice overview','Color',[1 1 1]);
end
