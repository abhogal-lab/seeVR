% Copyright (C) Alex A. Bhogal, 2025, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <GENERATEMOVIE: Montage movie of 4-D fMRI volumes + stimulus trace. >
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

% -------------------------------------------------------------------------
%   ▸ Displays each slice as     [ source | processed ]
%   ▸ Bottom row shows the stimulus trace with a moving time marker.
%   ▸ Encodes the montage as MP4 (default) or animated GIF.
%   ▸ A vertical colour-bar (%-signal-change) spans the right edge.
%
%───────────────────────────────────────────────────────────────────────────
% SYNTAX
%     generateMovie(mask,source,data,cmap,trace)           % all defaults
%     generateMovie(mask,source,data,cmap,trace,opts)      % customise
%
% REQUIRED INPUTS
%     mask      3-D logical [X×Y×Z]    – voxels to display
%     source    4-D numeric [X×Y×Z×T]  – raw volumes
%     data      4-D numeric [X×Y×Z×T]  – processed volumes (same geometry)
%     trace     1-D numeric (length=T) – stimulus / regressor / mean signal
%
% OPTIONAL STRUCT  OPTS            (all fields optional)
% ────────────────────────────────
%     type        0 = no file | 1 = GIF | 2 = MP4                 [2]
%     filename    output base name                                ['CVR_movie']
%     resultsdir  output folder                                   [pwd]
%     framerate   frames·s-1 for MP4                              [30]
%
%     Temporal controls
%     win         # dynamics averaged for processed view          [1]
%     start_ind   first dynamic rendered (1-based)                [1]
%
%     Slice controls
%     start       first slice index (1-based)                     [1]
%     step        slice increment                                 [1]
%     row         desired image rows   (≤6)                       [auto]
%     col         desired column-pairs (≤6)                       [auto]
%
%     Display controls
%     scale       [low high] for processed colour-limits          [-1 1]
%
% GRID LIMIT
%     Montage is capped at **6 image rows × 6 column-pairs**
%     (→ 12 subplot columns; capacity = 36 slices).  Extra slices are
%     ignored with a console warning.
%
% DEPENDENCIES
%     MATLAB R2019b+ (tiledlayout) • freezeColors (optional)
%     brewermap/cbrewer (optional – nice RdYlBu map; falls back to parula)
%
% AUTHORSHIP
%     Original  : Alex Bhogal   • a.bhogal@umcutrecht.nl
%
%__________________________________________________________________________

function generateMovie(mask,source,data,cmap,trace,opts)
%GENERATEMOVIE  Montage movie of 4-D fMRI volumes + stimulus trace.
%   Full documentation omitted here for brevity ─ nothing changed
%   except the colour-bar implementation.  See previous version for help.
%__________________________________________________________________________

%% 0 — defaults
% ----------------------------------------------------------------------
% Default-option helper
% Call this immediately after you receive the user-supplied OPTS struct.
% ----------------------------------------------------------------------

% ---------- file / output ---------------------------------------------
if ~isfield(opts,'type'),       opts.type       = 2;               end
if ~isfield(opts,'filename'),   opts.filename   = 'CVR_movie';     end
if ~isfield(opts,'resultsdir'), opts.resultsdir = pwd;            end
if ~isfield(opts,'framerate'),  opts.framerate  = 24;             end

% ---------- temporal controls -----------------------------------------
if ~isfield(opts,'win'),        opts.win        = 3;              end
if ~isfield(opts,'step'),       opts.step       = 3;       end
if ~isfield(opts,'start_ind'),  opts.start_ind  = 1;              end

% ---------- slice controls --------------------------------------------
if ~isfield(opts,'start'),      opts.start      = 10;              end
if ~isfield(opts,'slice_step'), opts.slice_step = 1;              end
if ~isfield(opts,'row'),        opts.row        = 4;             end
if ~isfield(opts,'col'),        opts.col        = 6;             end

% ---------- display ----------------------------------------------------
if ~isfield(opts,'scale'),      opts.scale      = [-5 5];         end
win  = opts.win;
dyn  = size(data,4);

%% 1 — rotate & mask
mask      = logical(imrotate(mask,90));
source    = imrotate(source,90);   source(source==0)=NaN;
data      = imrotate(data,90);     data(data==0)=NaN;
if isempty(cmap); cmap = flip(brewermap(128,'spectral')); end
%% 2 — grid size  (≤6 rows × 6 pairs)
maxRows=6; maxPairs=6;
sliceIdx = opts.start:opts.step:size(source,3);
nslices  = numel(sliceIdx);

row_img  = ~isempty(opts.row)*min(maxRows,opts.row) ...
    + isempty(opts.row)*min(maxRows,ceil(nslices/maxPairs));
colPairs = ~isempty(opts.col)*min(maxPairs,opts.col) ...
    + isempty(opts.col)*min(maxPairs,ceil(nslices/row_img));

capacity = row_img*colPairs;
if nslices>capacity
    warning('generateMovie:GridCapacity',...
        'Showing only first %d of %d slices (grid limit %d×%d pairs).',...
        capacity,nslices,maxRows,maxPairs);
    sliceIdx = sliceIdx(1:capacity); nslices=capacity;
end
opts.row=row_img; opts.col=colPairs*2; traceRow = opts.row+1;

%% 3 — output initialisation  (safe close on exit)
v=[]; cleanupVideo=onCleanup(@()safeClose(v));
switch opts.type
    case 2
        v = VideoWriter(fullfile(opts.resultsdir,[opts.filename '.mp4']),'MPEG-4');
        v.FrameRate = opts.framerate; open(v)
    case 1
        gifFile = fullfile(opts.resultsdir,[opts.filename '.gif']);
end

%% 4 — main loop
fig=figure('Color','k','Units','pixels','Position',[10 80 1600 900]);

for ii=opts.start_ind:opts.step_size:dyn-win
    t = tiledlayout(traceRow,opts.col,'Padding','loose','TileSpacing','compact');
    sliceCounter=1;

    % 4a — images --------------------------------------------------------
    for r=1:opts.row
        for p=1:colPairs
            if sliceCounter>nslices, break, end
            s=sliceIdx(sliceCounter);

            nexttile;                    % source
            imagesc(mask(:,:,s).*source(:,:,s,ii));
            colormap(gray); if exist('freezeColors','file'), freezeColors,end
            axis off

            nexttile;                    % processed
            tmp = mean(data(:,:,s,ii:ii+win),4,'omitnan').*mask(:,:,s);
            imagesc(tmp,opts.scale); colormap(cmap);
            if exist('freezeColors','file'), freezeColors,end
            axis off

            sliceCounter=sliceCounter+1;
        end
    end

    % 4b — custom-positioned trace panel with colormap-shaded background ----
    % Get position of full grid to anchor custom axis
    parentPos = get(t, 'Position');         % [x y width height] of tiledlayout
    traceHeightFrac = 0.12;                 % adjustable height for trace subplot

    % Create trace axis manually (not part of tiledlayout)
    traceAx = axes('Parent', fig, ...
        'Position', [parentPos(1), ...
        parentPos(2), ...
        parentPos(3), ...
        parentPos(4) * traceHeightFrac], ...
        'Color', 'k');

    hold(traceAx, 'on')

    % Normalize trace to [0,1] for shading
    traceNorm = (trace - min(trace)) / (max(trace) - min(trace) + eps);

    % Define background colormap for the trace panel
    if isfield(opts, 'trace_cmap') && ~isempty(opts.trace_cmap)
        traceCmap = opts.trace_cmap;  % User-supplied colormap (Nx3)
    else
        traceCmap = gray(128);        % Default to grayscale if not specified
    end

    % Interpolate trace values to colormap indices
    cmapIdx = round(1 + traceNorm * (size(traceCmap, 1)-1));
    cmapIdx = min(max(cmapIdx, 1), size(traceCmap,1));

    % Draw shaded background rectangles
    for tIdx = 1:dyn
        fill(traceAx, ...
            [tIdx-0.5, tIdx+0.5, tIdx+0.5, tIdx-0.5], ...
            [0, 0, 1, 1], ...
            traceCmap(cmapIdx(tIdx),:), ...
            'EdgeColor', 'none');
    end

    % Plot trace line (top color of main cmap)
    plot(traceNorm, 'Parent', traceAx, ...
        'LineWidth', 2.5, 'Color', cmap(end,:));

    % Vertical time marker (bottom color of main cmap)
    tMid = ii + (win-1)/2;
    line(traceAx, [tMid tMid], [0 1], ...
        'Color', cmap(1,:), 'LineWidth', 3, 'Clipping', 'off');

    % Axis formatting
    xlim(traceAx, [1 dyn]);
    ylim(traceAx, [0 1]);
    set(traceAx, 'Color', 'k', 'XColor', 'w', 'YColor', 'none', ...
        'LineWidth', 1.5, 'FontSize', 10, ...
        'XTick', [], 'YTick', []);
    xlabel(traceAx, 'Dynamic scans', 'Color', 'w');


    % 4c — colour-bar (right edge) -------------------------------------
    % Attach to a dummy axes for controlled positioning
    cbAx = axes('Parent',fig,'Position',[0.93 0.15 0.02 0.7],'Visible','off');
    colormap(cbAx, cmap);
    caxis(cbAx, opts.scale);
    cb = colorbar(cbAx);
    cb.Units = 'normalized';
    cb.Position = [0.95 0.15 0.015 0.7];  % [x y width height]
    cb.Label.String = '% signal change';
    cb.Color = 'w';
    cb.FontSize = 12;

    % 4d — write frame ---------------------------------------------------
    drawnow nocallbacks
    frame=getframe(fig);

    switch opts.type
        case 2, writeVideo(v,frame);
        case 1
            [imind,cm]=rgb2ind(frame2im(frame),256);
            if ii==opts.start_ind
                imwrite(imind,cm,gifFile,'gif','Loopcount',inf,'DelayTime',0);
            else
                imwrite(imind,cm,gifFile,'gif','WriteMode','append','DelayTime',0);
            end
    end
    if exist('cb','var'), delete(cb); end
    if exist('cbAx','var'), delete(cbAx); end
    delete(t);  % safe to delete the tiledlayout always
end

if opts.type==2, close(v), end
close(fig)
end

%% helpers ---------------------------------------------------------------
function s=copyfields(s,dflt)
fn=fieldnames(dflt);
for k=1:numel(fn)
    if ~isfield(s,fn{k})||isempty(s.(fn{k})), s.(fn{k})=dflt.(fn{k}); end
end
end
function safeClose(w)
if ~isempty(w)&&isvalid(w)&&strcmp(w.Status,'open'), try close(w), end, end
end
