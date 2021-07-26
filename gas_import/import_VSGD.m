%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [MRTimes2, G1targetvolumemL, G1FCO1, G1FO1, G2FCO1, G2FO1, G1measuredbreathmL, G2measuredbreathmL, BlendermeasuredtotalmL, BlendermeasuredCO2mL, BlendermeasuredO2mL, BlendermeasuredS1mL, BlendermeasuredS2mL] = import_VSGD(filename, dataLines)

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["MRTimes2", "G1targetvolumemL", "G1FCO1", "G1FO1", "G2FCO1", "G2FO1", "G1measuredbreathmL", "G2measuredbreathmL", "BlendermeasuredtotalmL", "BlendermeasuredCO2mL", "BlendermeasuredO2mL", "BlendermeasuredS1mL", "BlendermeasuredS2mL"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
MRTimes2 = tbl.MRTimes2;
G1targetvolumemL = tbl.G1targetvolumemL;
G1FCO1 = tbl.G1FCO1;
G1FO1 = tbl.G1FO1;
G2FCO1 = tbl.G2FCO1;
G2FO1 = tbl.G2FO1;
G1measuredbreathmL = tbl.G1measuredbreathmL;
G2measuredbreathmL = tbl.G2measuredbreathmL;
BlendermeasuredtotalmL = tbl.BlendermeasuredtotalmL;
BlendermeasuredCO2mL = tbl.BlendermeasuredCO2mL;
BlendermeasuredO2mL = tbl.BlendermeasuredO2mL;
BlendermeasuredS1mL = tbl.BlendermeasuredS1mL;
BlendermeasuredS2mL = tbl.BlendermeasuredS2mL;
end