% ===========================================================
% File: createParallelProgressBar.m
% Description: Parfor progress bar using DataQueue + waitbar
%              (NO Java/JavaPeer usage). Supports chunked updates.
% License: MIT
% ===========================================================

function queue = createParallelProgressBar(totalIterations, varargin)
% createParallelProgressBar  Progress bar for parfor using DataQueue (no Java).
%
% queue = createParallelProgressBar(totalIterations)
% queue = createParallelProgressBar(totalIterations,'Message','Processing...','Name','Computation Progress')
% queue = createParallelProgressBar(totalIterations,'UpdateEvery',25,'ColorStart',[...],'ColorEnd',[...])
%
% In parfor (recommended for huge N):
%   chunk = 250;
%   local = 0;
%   parfor i = 1:N
%       % ... work ...
%       local = local + 1;
%       if local >= chunk
%           send(queue, local);   % increment by local (chunk)
%           local = 0;
%       end
%   end
%   % (Optionally flush remainder per worker; see notes below.)
%
% Notes:
% - Uses parallel.pool.DataQueue + afterEach (client-side UI updates).
% - Uses MATLAB HG objects to color the waitbar (no JavaPeer).
% - Each send(queue, n) increments progress by n (default n=1).
% - If the user closes the waitbar, updates stop safely.
%
% Optional: For exact completion with chunking, flush remainder per worker
% using an spmd block after the parfor (each worker holds its own "local"
% counter; you can implement with persistent per-worker counters if needed).

    % -------- Parse options --------
    ip = inputParser;
    ip.addRequired('totalIterations', @(x) validateattributes(x, {'numeric'}, {'scalar','positive','integer'}));
    ip.addParameter('Message', 'Processing...', @(s) ischar(s) || isstring(s));
    ip.addParameter('Name', 'Computation Progress', @(s) ischar(s) || isstring(s));
    ip.addParameter('ColorStart', [171, 94, 0]/255, @(v) isnumeric(v) && numel(v)==3);
    ip.addParameter('ColorEnd',   [12, 123, 220]/255, @(v) isnumeric(v) && numel(v)==3);
    ip.addParameter('UpdateEvery', 1, @(x) validateattributes(x, {'numeric'}, {'scalar','positive','integer'})); % throttle UI redraws
    ip.parse(totalIterations, varargin{:});
    opts = ip.Results;

    % -------- Create queue + waitbar (client side) --------
    queue = parallel.pool.DataQueue;

    h = waitbar(0, char(opts.Message), 'Name', char(opts.Name));
    drawnow;

    % Try to find the filled bar patch for coloring
    barPatch = localFindWaitbarPatch(h);
    if ~isempty(barPatch) && isgraphics(barPatch)
        try
            barPatch.FaceColor = opts.ColorStart;
            barPatch.EdgeColor = 'none';
        catch
            % Ignore styling failures (varies by MATLAB release)
        end
    end

    % Store mutable state in appdata
    state = struct();
    state.Total       = double(totalIterations);
    state.Count       = 0;
    state.UpdateEvery = double(opts.UpdateEvery);
    state.UserClosed  = false;
    state.LastDrawn   = 0; % last count when UI was updated
    setappdata(h, 'ParallelProgressBarState', state);

    % If user closes figure early, stop updating safely
    h.CloseRequestFcn = @(src,evt) localUserClose(src);

    % Listener must be kept alive: store in appdata to avoid GC
    L = afterEach(queue, @(n) localUpdate(h, barPatch, opts, n));
    setappdata(h, 'ParallelProgressBarListener', L);
end

% ===== Client-side update callback =====
function localUpdate(h, barPatch, opts, n)
    if ~isgraphics(h)
        return; % already closed
    end

    state = getappdata(h, 'ParallelProgressBarState');
    if isempty(state) || state.UserClosed
        if isgraphics(h), delete(h); end
        return;
    end

    if nargin < 4 || isempty(n)
        n = 1;
    end
    n = double(n);

    % Increment count by payload (supports chunking)
    state.Count = state.Count + n;
    if state.Count > state.Total
        state.Count = state.Total;
    end

    % Throttle UI updates (redraw every UpdateEvery increments)
    if (state.Count - state.LastDrawn) < state.UpdateEvery && state.Count < state.Total
        setappdata(h, 'ParallelProgressBarState', state);
        return;
    end
    state.LastDrawn = state.Count;

    frac = min(1, state.Count / state.Total);
    pct  = 100 * frac;

    msg = sprintf('%s  (%d/%d, %5.1f%%)', char(opts.Message), ...
        round(state.Count), round(state.Total), pct);

    % Update waitbar
    try
        waitbar(frac, h, msg);
    catch
        return; % figure may have closed
    end

    % Update bar color (HG patch)
    if ~isempty(barPatch) && isgraphics(barPatch)
        c = (1 - frac) * opts.ColorStart + frac * opts.ColorEnd;
        try
            barPatch.FaceColor = c;
        catch
            % ignore (release differences)
        end
    end

    drawnow limitrate nocallbacks;

    % Auto-close on completion
    if state.Count >= state.Total
        if isgraphics(h), delete(h); end
    else
        setappdata(h, 'ParallelProgressBarState', state);
    end
end

% ===== Find waitbar filled-bar patch (graphics, not Java) =====
function barPatch = localFindWaitbarPatch(h)
    barPatch = [];
    if ~isgraphics(h), return; end

    patches = findall(h, 'Type', 'patch');
    if isempty(patches), return; end

    % Heuristic: pick a patch that looks like a rectangle
    for k = 1:numel(patches)
        pk = patches(k);
        try
            xd = pk.XData; yd = pk.YData;
            if isnumeric(xd) && isnumeric(yd) && numel(xd) >= 4 && numel(yd) >= 4
                barPatch = pk;
                return;
            end
        catch
        end
    end

    % Fallback
    barPatch = patches(1);
end

% ===== Handle user closing the waitbar =====
function localUserClose(h)
    if ~isgraphics(h), return; end
    state = getappdata(h, 'ParallelProgressBarState');
    if ~isempty(state)
        state.UserClosed = true;
        setappdata(h, 'ParallelProgressBarState', state);
    end
    delete(h);
end
