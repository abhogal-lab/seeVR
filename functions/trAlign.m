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
handles.varCheck = length(varargin);
handles.output = hObject;
handles.corrvec = varargin{1,1}; if iscolumn(handles.corrvec); handles.corrvec = handles.corrvec'; end
handles.scale = mean(handles.corrvec);

if handles.varCheck > 3
handles.corrvec2 = varargin{1,2}; if iscolumn(handles.corrvec2); handles.corrvec2 = handles.corrvec2'; end
handles.TS = varargin{1,3}; if iscolumn(handles.TS); handles.TS = handles.TS'; end
else
handles.TS = varargin{1,2}; if iscolumn(handles.TS); handles.TS = handles.TS'; end
end

handles.opts = varargin{1,end};

if isfield(handles.opts, 'figdir')
handles.figdir = handles.opts.figdir;
elseif isfield(handles.opts, 'savedir')
handles.figdir = handles.opts.savedir; 
else
handles.figdir = [pwd,'/'];
end

handles.dyn = length(handles.TS);
[handles.r,handles.lags] = xcorr(rescale(handles.corrvec),rescale(handles.TS));
handles.lagmask = handles.lags; handles.lagmask(handles.lags < 0) = 0; handles.lagmask(handles.lagmask > 0) = 1; 
handles.r = handles.r.*handles.lagmask;
[handles.M,handles.I] = max((handles.r(:)));
handles.maxcorr = handles.lags(handles.I);
handles.xdata = 1:1:handles.dyn;

handles.midlag = find(handles.lags == 0);
handles.st = handles.maxcorr;
handles.en = handles.st + handles.dyn -1;
handles.probe1 = handles.corrvec(1,handles.st:handles.en);

if handles.varCheck > 3
handles.probe2 = handles.corrvec2(1,handles.st:handles.en);
end

%checkplots
%figure; plot(rescale(handles.probe1)); hold on; plot(rescale(handles.TS)); if handles.varCheck > 3; plot(rescale(handles.probe2)); end 

%setup slider to allow translation

minRange = handles.maxcorr; 
maxRange = length(handles.corrvec) - length(handles.TS) - handles.maxcorr;
slRange = maxRange + minRange;
corrWin = 2*maxRange;

handles.lastSliderVal = 0; handles.offset = handles.maxcorr;

set(handles.slider1, 'Min', -(minRange-1), 'Max', maxRange , 'Sliderstep', [1/slRange 1/slRange] ,'Value', handles.lastSliderVal);

% if (corrWin + handles.midlag) > length(handles.lags)
% handles.window = [handles.midlag+1:1:length(handles.lags)];
% else
% handles.window = [handles.midlag+1:1:(handles.midlag + corrWin)];
% end
%plot initial correlation
handles.stPt = handles.maxcorr;
handles.window =  [handles.midlag+1:1:length(handles.lags)];


axes(handles.traces)
[ax,h1,h2] = plotyy(handles.xdata,handles.probe1,handles.xdata,handles.TS);
set(ax,'xtick',[]); set(ax,'xticklabel',[]); set(ax,'ytick',[]); set(ax,'yticklabel',[]); set(ax,'xlim',[0 handles.dyn]);
set(h1,'color','m','Linewidth', 1.5); set(h2,'color','k','Linewidth', 1.5); set(ax,'xlim',[0 handles.dyn]);
axes(handles.correlation)
plot(handles.lags(1,handles.window),handles.r(1,handles.window),'-k', 'Linewidth', 1.5);
set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[])
hold on; xline(handles.stPt,'-r', 'Linewidth', 2);

% Update handles structure

guidata(hObject, handles);

% UIWAIT makes trAlign wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trAlign_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% Get default command line output from handles structure
%varargout{1} = handles.output;



% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)

handles.probe1 = [];
handles.probe2 = [];
handles.lastSliderVal = round(get(hObject,'Value'));
handles.offset = round(get(hObject,'Value'));
handles.st = handles.maxcorr + handles.offset;
handles.en = handles.st + handles.dyn;
handles.probe1 = handles.corrvec(1,handles.st:handles.en-1);
if handles.varCheck > 3
handles.probe2 = handles.corrvec2(1,handles.st:handles.en-1);
end
disp(['offset is: ',int2str(handles.offset)])

if (handles.stPt + handles.offset) > length(handles.window)
handles.window = [handles.midlag+1:1:(handles.midlag + (handles.stPt + ceil(handles.dyn/2)))];
end


axes(handles.traces);
[ax,h1,h2] = plotyy(handles.xdata,handles.probe1,handles.xdata,handles.TS);
set(ax,'xtick',[]); set(ax,'xticklabel',[]); set(ax,'ytick',[]); set(ax,'yticklabel',[]); set(ax,'xlim',[0 handles.dyn]);
set(h1,'color','m','Linewidth', 1.5); set(h2,'color','k','Linewidth', 1.5); set(ax,'xlim',[0 handles.dyn]);
axes(handles.correlation);
cla()
plot(handles.lags(1,handles.window),handles.r(1,handles.window),'-k', 'Linewidth', 2);
hold on; xline(handles.stPt + handles.offset,'-r', 'Linewidth', 1.5); hold off;
set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[])
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

%save traces
assignin('base', 'offset',  handles.offset)
assignin('base', 'probe1',  handles.probe1)

if handles.varCheck > 3
assignin('base', 'probe2',  handles.probe2)
nrplots = 3;
else
nrplots = 1;
end

xdata = [1:1:length(handles.TS)];
%save([handles.figdir,'traces.mat'], 'handles.probe1', 'handles.probe2')
figure(11); 
if nrplots == 3
subplot(nrplots,1,1); [ax,h1,h2] = plotyy(xdata,handles.probe1,xdata,handles.TS);
set(h1,'color','m'); set(h2,'color','k'); set(ax,'xlim',[0 handles.dyn]);
set(ax(1),'ycolor','m'); set(ax(2),'ycolor','k');
xlabel('scan number'); ylabel('probe'); title('main probe & signal');
subplot(nrplots,1,2) ; [ax2,h1,h2] = plotyy(xdata,handles.probe2,xdata,handles.TS);
set(h1,'color','b'); set(h2,'color','k'); set(ax2,'xlim',[0 handles.dyn]);
xlabel('scan number'); ylabel('probe2'); title('secondary probe & signal');
set(ax2(1),'ycolor','b'); set(ax2(2),'ycolor','k');
subplot(nrplots,1,3); [ax3,h1,h2] = plotyy(xdata,handles.probe2,xdata,handles.probe1);
set(h1,'color','b'); set(h2,'color','m'); set(ax3,'xlim',[0 handles.dyn]);
set(ax3(1),'ycolor','b'); set(ax3(2),'ycolor','m');
title('probes');
else
[ax,h1,h2] = plotyy(xdata,handles.probe1,xdata,handles.TS);
set(h1,'color','m'); set(h2,'color','k'); set(ax,'xlim',[0 handles.dyn]);
set(ax(1),'ycolor','m'); set(ax(2),'ycolor','k');
xlabel('scan number'); ylabel('probe'); title('main probe & signal');    
end
    
    
saveas(gcf,fullfile(handles.figdir,'trace_alignment.fig'));
disp('Closing alignment tool to continue analysis')
uiresume(handles.figure1);
guidata(hObject, handles);
close trAlign
