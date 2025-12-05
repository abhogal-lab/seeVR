% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <trAlign: aligns an input signal with a reference signal >
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
%
%
% inputs: this function takes either 3 or 4 inputs. The first input is the
% trace to be aligned. This is typically the end-tidal CO2 trace that has
% been resampled to the TR of the scan (use resampletoTR.m). A second
% physiolgical trace can be passed as the second input (for example petO2)
% to be aligned as a 'passenger trace'. The third input is the reference
% signal. For example the meanTimeseries in the GMmask. If only 3 inputs
% are supplied, then the first input should be the physiological trace and
% the second input the reference signal. The last input should be the opts
% structure.
%
% usage: trAlign(petCO2, petO2, reference, opts)
%        trAlign(petCO2, reference, opts)
%
% Outputs (optional, in addition to base workspace variables):
%   [probe1, probe2, probe3, probe4] = trAlign(...)
%   probe2 / probe4 are [] if no passenger trace is provided.

function varargout = trAlign(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @trAlign_OpeningFcn, ...
    'gui_OutputFcn',  @trAlign_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before trAlign is made visible.
function trAlign_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.

% Parse inputs more explicitly
nIn = numel(varargin);
handles.varCheck = nIn;

if nIn == 3
    % trAlign(probe1, TS, opts)
    handles.corrvec = varargin{1};
    handles.TS      = varargin{2};
    handles.opts    = varargin{3};
elseif nIn == 4
    % trAlign(probe1, probe2, TS, opts)
    handles.corrvec  = varargin{1};
    handles.corrvec2 = varargin{2};
    handles.TS       = varargin{3};
    handles.opts     = varargin{4};
else
    error('trAlign:WrongNArgs', ...
        'Usage: trAlign(probe1, TS, opts) or trAlign(probe1, probe2, TS, opts)');
end

% Ensure row vectors
if iscolumn(handles.corrvec),  handles.corrvec  = handles.corrvec.';  end
if iscolumn(handles.TS),       handles.TS       = handles.TS.';       end
if isfield(handles,'corrvec2') && iscolumn(handles.corrvec2)
    handles.corrvec2 = handles.corrvec2.';
end

handles.output = hObject;

% Figure directory
if isfield(handles.opts, 'figdir')
    handles.figdir = handles.opts.figdir;
elseif isfield(handles.opts, 'savedir')
    handles.figdir = handles.opts.savedir;
else
    handles.figdir = pwd;
end

handles.dyn   = numel(handles.TS);       % TS length is ground truth
handles.xdata = 1:handles.dyn;

% Pre-smooth master signals (we will always slice from these for display)
handles.corrvec_sm = smooth(handles.corrvec(:)).';
handles.TS_sm      = smooth(handles.TS(:)).';

if isfield(handles,'corrvec2')
    handles.corrvec2_sm = smooth(handles.corrvec2(:)).';
end

% Cross-correlation (full lags, for plotting + negative lags)
[handles.rAll, handles.lags] = xcorr(rescale(handles.corrvec), ...
                                     rescale(handles.TS), ...
                                     'none');

% ---------- choose initial lag: best that allows a full window ----------
posIdx  = handles.lags >= 0;
rPos    = handles.rAll(posIdx);
lagsPos = handles.lags(posIdx);

NcorrSm = numel(handles.corrvec_sm);
dyn     = handles.dyn;

% Sort positive-lag correlations by magnitude (descending)
[sortedR, sortIdx] = sort(rPos, 'descend');
sortedLags         = lagsPos(sortIdx);

chosenLag = sortedLags(1);
chosenR   = sortedR(1);

for k = 1:numel(sortedLags)
    L  = sortedLags(k);
    st = L;               % original semantics: st = lag
    en = st + dyn - 1;
    if st >= 1 && en <= NcorrSm
        chosenLag = L;
        chosenR   = sortedR(k);
        break;
    end
end

handles.maxcorr = chosenLag;       % chosen initial lag
handles.M       = chosenR;
handles.midlag  = find(handles.lags == 0, 1);  % index of zero lag

% Store r for plotting (now includes both negative and positive lags)
handles.r = handles.rAll;

% Initial window for corr plot: show full lag range
handles.window = 1:numel(handles.lags);

% Initial probe window based on maxcorr (keep original semantics)
handles.st = handles.maxcorr;
handles.en = handles.st + handles.dyn - 1;

Ncorr = numel(handles.corrvec_sm);

% Ensure that we can get a window of length dyn from corrvec_sm if possible
if Ncorr >= handles.dyn
    if handles.st < 1
        handles.st = 1;
    end
    if handles.st > Ncorr - handles.dyn + 1
        handles.st = Ncorr - handles.dyn + 1;
    end
    handles.en = handles.st + handles.dyn - 1;
else
    % corrvec shorter than TS: take full corrvec and pad probe to match dyn
    handles.st = 1;
    handles.en = Ncorr;
end

% Main aligned probes (smoothed for display) - always pad to dyn
win1 = handles.corrvec_sm(handles.st:handles.en);
if numel(win1) < handles.dyn
    win1 = [win1, nan(1, handles.dyn - numel(win1))];
end
handles.probe1 = win1(1:handles.dyn);

if isfield(handles,'corrvec2_sm')
    N2  = numel(handles.corrvec2_sm);
    idx = handles.st:handles.en;
    vals = nan(1, numel(idx));
    valid = idx >= 1 & idx <= N2;
    vals(valid) = handles.corrvec2_sm(idx(valid));
    if numel(vals) < handles.dyn
        vals = [vals, nan(1, handles.dyn - numel(vals))];
    end
    handles.probe2 = vals(1:handles.dyn);
end

% Extended probes (±5% extra from original input corrvec/corrvec2)
handles = computeExtendedProbes(handles);

% ------------------------ Slider setup ------------------------
% offset 0 = auto-lag guess, slider refines around it.
leftMaxShift  = max(0, handles.maxcorr - 1);
rightMaxShift = max(0, numel(handles.corrvec) - numel(handles.TS) - handles.maxcorr);

sliderMin = -leftMaxShift;
sliderMax =  rightMaxShift;

% If range is degenerate (no real room to move), still give a tiny range
if sliderMax <= sliderMin
    sliderMin = -1;
    sliderMax =  1;
end

slRange = sliderMax - sliderMin;
if slRange <= 0
    slRange = 1;
end

handles.lastSliderVal = 0;
handles.offset        = 0;
handles.stPt          = handles.maxcorr;  % initial lag indicator for plot

set(handles.slider1, ...
    'Min', sliderMin, ...
    'Max', sliderMax, ...
    'Sliderstep', [1/slRange 1/slRange], ...
    'Value', handles.lastSliderVal);

% --- Initial plotting ---

% Traces panel: probe1 vs TS
axes(handles.traces); cla(handles.traces);
[ax,h1,h2] = plotyy(handles.xdata, handles.probe1, ...
                    handles.xdata, handles.TS_sm);
set(ax,'XTick',[], 'XTickLabel',[], 'YTick',[], 'YTickLabel',[]);
set(ax, 'XLim', [0, handles.dyn]);
set(h1, 'Color','m', 'LineWidth', 1.5, 'HitTest','off');
set(h2, 'Color','k', 'LineWidth', 1.5, 'HitTest','off');
xlabel(handles.traces,'Click to place reference line','FontSize',8);

% Allow clicking in traces axes to add a vertical reference line
set(handles.traces, 'ButtonDownFcn', @traces_ButtonDownFcn);

% Correlation panel: full lag range, with initial peak lag and current lag
axes(handles.correlation); cla(handles.correlation);
plot(handles.lags(handles.window), handles.r(handles.window), '-k', 'LineWidth', 1.5);
hold on;

curLag = handles.stPt + handles.offset;  % at opening, offset=0 -> maxcorr
[curR, found] = getCorrAtLag(handles.rAll, handles.lags, curLag);
if ~found
    curR = NaN;
end

% Vertical line at initial correlation lag (maxcorr) - blue dashed
xline(handles.stPt, '--b', 'LineWidth', 1.5);

% Vertical line at current lag (red)
xline(curLag, '-r', 'LineWidth', 2);

set(gca,'XTick',[], 'XTickLabel',[], 'YTick',[], 'YTickLabel',[]);
title('xcorr', 'FontWeight','normal');

% Move numeric info to x-label area: between plots and slider
xlabel(sprintf('lag = %d   |   offset = %d   |   xcorr = %.3f', ...
               curLag, handles.offset, curR));

hold off;

% Save and wait for user interaction
guidata(hObject, handles);
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trAlign_OutputFcn(hObject, eventdata, handles)
% Now optionally return aligned and extended probes, in addition
% to still exporting to base in pushbutton1_Callback.

if isempty(handles)
    % User closed manually
    for k = 1:nargout
        varargout{k} = [];
    end
    return;
end

% Default first output is still the GUI handle (for backwards compat),
% but if the user asks for multiple outputs, give probes instead.
if nargout == 1
    varargout{1} = handles.output;
elseif nargout >= 1
    % 1: probe1
    varargout{1} = handles.probe1;
    % 2: probe2 (if present)
    if nargout >= 2
        if isfield(handles,'probe2')
            varargout{2} = handles.probe2;
        else
            varargout{2} = [];
        end
    end
    % 3: probe3 (extended main)
    if nargout >= 3
        if isfield(handles,'probe3')
            varargout{3} = handles.probe3;
        else
            varargout{3} = [];
        end
    end
    % 4: probe4 (extended passenger)
    if nargout >= 4
        if isfield(handles,'probe4')
            varargout{4} = handles.probe4;
        else
            varargout{4} = [];
        end
    end
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)

handles.lastSliderVal = round(get(hObject,'Value'));
handles.offset        = handles.lastSliderVal;

% Compute new start/end indices for probe window
handles.st = handles.maxcorr + handles.offset;
handles.en = handles.st + handles.dyn - 1;

Ncorr = numel(handles.corrvec_sm);

if Ncorr >= handles.dyn
    if handles.st < 1
        handles.st = 1;
    end
    if handles.st > Ncorr - handles.dyn + 1
        handles.st = Ncorr - handles.dyn + 1;
    end
    handles.en = handles.st + handles.dyn - 1;
else
    handles.st = 1;
    handles.en = Ncorr;
end

% Update main probes (smoothed, padded to dyn)
win1 = handles.corrvec_sm(handles.st:handles.en);
if numel(win1) < handles.dyn
    win1 = [win1, nan(1, handles.dyn - numel(win1))];
end
handles.probe1 = win1(1:handles.dyn);

if isfield(handles,'corrvec2_sm')
    N2  = numel(handles.corrvec2_sm);
    idx = handles.st:handles.en;
    vals = nan(1, numel(idx));
    valid = idx >= 1 & idx <= N2;
    vals(valid) = handles.corrvec2_sm(idx(valid));
    if numel(vals) < handles.dyn
        vals = [vals, nan(1, handles.dyn - numel(vals))];
    end
    handles.probe2 = vals(1:handles.dyn);
end

% Update extended probes (from original input with ±5% extra)
handles = computeExtendedProbes(handles);

disp(['offset is: ', int2str(handles.offset)]);

% --- Update traces plot ---
axes(handles.traces); cla(handles.traces);
[ax,h1,h2] = plotyy(handles.xdata, handles.probe1, ...
                    handles.xdata, handles.TS_sm);
set(ax,'XTick',[], 'XTickLabel',[], 'YTick',[], 'YTickLabel',[]);
set(ax, 'XLim', [0, handles.dyn]);
set(h1, 'Color','m', 'LineWidth', 1.5, 'HitTest','off');
set(h2, 'Color','k', 'LineWidth', 1.5, 'HitTest','off');
xlabel(handles.traces,'Click to place reference line','FontSize',8);

% Re-draw reference line in traces plot if it exists
if isfield(handles,'refX')
    hold on;
    handles.refLine = xline(handles.refX, '--g', 'LineWidth', 1.5);
    set(handles.refLine,'HitTest','off');
    hold off;
end

% Make sure clicking in traces axes still works
set(handles.traces, 'ButtonDownFcn', @traces_ButtonDownFcn);

% --- Update correlation plot with current lag + r value + offset ---
axes(handles.correlation); cla(handles.correlation);
plot(handles.lags(handles.window), handles.r(handles.window), '-k', 'LineWidth', 2);
hold on;

curLag = handles.stPt + handles.offset;
[curR, found] = getCorrAtLag(handles.rAll, handles.lags, curLag);
if ~found
    curR = NaN;
end

% Initial lag (blue dashed)
xline(handles.stPt, '--b', 'LineWidth', 1.5);

% Current lag (red)
xline(curLag, '-r', 'LineWidth', 1.5);

set(gca,'XTick',[], 'XTickLabel',[], 'YTick',[], 'YTickLabel',[]);
title('xcorr', 'FontWeight','normal');
xlabel(sprintf('lag = %d   |   offset = %d   |   xcorr = %.3f', ...
               curLag, handles.offset, curR));

hold off;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), ...
           get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

% Export aligned traces & offset to base workspace (keep original behaviour)
assignin('base', 'offset',  handles.offset);
assignin('base', 'probe1',  handles.probe1);

if handles.varCheck > 3 && isfield(handles,'probe2')
    assignin('base', 'probe2', handles.probe2);
end

% Also export extended probes (±5% extra)
if isfield(handles,'probe3')
    assignin('base', 'probe3', handles.probe3);
end
if isfield(handles,'probe4')
    assignin('base', 'probe4', handles.probe4);
end

% Decide number of summary plots based on presence of passenger probe
if handles.varCheck > 3 && isfield(handles,'probe2')
    nrplots = 3;
else
    nrplots = 1;
end

xdata = 1:numel(handles.TS_sm);

% Diagnostic figure with alignment summary
figure(11); clf;

if nrplots == 3
    % 1) main probe vs TS
    subplot(nrplots,1,1);
    [ax,h1,h2] = plotyy(xdata, handles.probe1, xdata, handles.TS_sm);
    set(h1,'Color','m'); set(h2,'Color','k');
    set(ax,'XLim',[0 handles.dyn]);
    set(ax(1),'YColor','m'); set(ax(2),'YColor','k');
    xlabel('scan number'); ylabel('probe');
    title('Main probe & signal');

    % 2) secondary probe vs TS
    subplot(nrplots,1,2);
    [ax2,h3,h4] = plotyy(xdata, handles.probe2, xdata, handles.TS_sm);
    set(h3,'Color','b'); set(h4,'Color','k');
    set(ax2,'XLim',[0 handles.dyn]);
    set(ax2(1),'YColor','b'); set(ax2(2),'YColor','k');
    xlabel('scan number'); ylabel('probe2');
    title('Secondary probe & signal');

    % 3) probes vs each other
    subplot(nrplots,1,3);
    [ax3,h5,h6] = plotyy(xdata, handles.probe2, xdata, handles.probe1);
    set(h5,'Color','b'); set(h6,'Color','m');
    set(ax3,'XLim',[0 handles.dyn]);
    set(ax3(1),'YColor','b'); set(ax3(2),'YColor','m');
    xlabel('scan number'); ylabel('probes');
    title('Probes (main vs secondary)');

else
    % Only main probe vs TS
    [ax,h1,h2] = plotyy(xdata, handles.probe1, xdata, handles.TS_sm);
    set(h1,'Color','m'); set(h2,'Color','k');
    set(ax,'XLim',[0 handles.dyn]);
    set(ax(1),'YColor','m'); set(ax(2),'YColor','k');
    xlabel('scan number'); ylabel('probe');
    title('Main probe & signal');
end

% Save figure
saveas(gcf, fullfile(handles.figdir,'trace_alignment.fig'));

disp('Closing alignment tool to continue analysis');
uiresume(handles.figure1);
guidata(hObject, handles);
close trAlign;


% === Helper: get correlation value at a given lag ========================
function [val, found] = getCorrAtLag(rAll, lags, lagQuery)
% Return correlation value at exact lag if present

idx = find(lags == lagQuery, 1);
if isempty(idx)
    val   = NaN;
    found = false;
else
    val   = rAll(idx);
    found = true;
end


% === Helper: compute extended probes (±5% of dyn) =======================
function handles = computeExtendedProbes(handles)

extra = max(1, round(0.05 * handles.dyn));

N1 = numel(handles.corrvec);
stEx = max(handles.st - extra, 1);
enEx = min(handles.en + extra, N1);
handles.probe3 = handles.corrvec(stEx:enEx);   % raw main trace with margin

if isfield(handles,'corrvec2')
    N2 = numel(handles.corrvec2);
    % Use same relative window, but clip to passenger length
    stEx2 = max(handles.st - extra, 1);
    enEx2 = min(handles.en + extra, N2);
    handles.probe4 = handles.corrvec2(stEx2:enEx2);
end


% === Callback: mouse click on traces axes to add reference line =========
function traces_ButtonDownFcn(hAx, ~)

handles = guidata(hAx);
if isempty(handles)
    return;
end

% Get x-position of click in data coordinates (sample index)
cp = get(hAx, 'CurrentPoint');
xClick = round(cp(1,1));
xClick = max(1, min(handles.dyn, xClick));  % clamp to valid sample range
handles.refX = xClick;

% Redraw traces plot to ensure both signals stay visible
axes(handles.traces); cla(handles.traces);
[ax,h1,h2] = plotyy(handles.xdata, handles.probe1, ...
                    handles.xdata, handles.TS_sm);
set(ax,'XTick',[], 'XTickLabel',[], 'YTick',[], 'YTickLabel',[]);
set(ax, 'XLim', [0, handles.dyn]);
set(h1, 'Color','m', 'LineWidth', 1.5, 'HitTest','off');
set(h2, 'Color','k', 'LineWidth', 1.5, 'HitTest','off');
xlabel(handles.traces,'Click to place reference line','FontSize',8);

hold on;
handles.refLine = xline(handles.refX, '--g', 'LineWidth', 1.5);
set(handles.refLine,'HitTest','off');
hold off;

% Make sure axes keeps the click callback
set(handles.traces, 'ButtonDownFcn', @traces_ButtonDownFcn);

guidata(hAx, handles);
