%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [MRTimes3, CtrlRoomTimes, Event] = import_Events(filename, dataLines)

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["MRTimes3", "CtrlRoomTimes", "Event"];
opts.VariableTypes = ["double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Event", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Event", "EmptyFieldRule", "auto");

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
MRTimes3 = tbl.MRTimes3;
CtrlRoomTimes = tbl.CtrlRoomTimes;
Event = tbl.Event;
end