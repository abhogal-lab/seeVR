% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <import_VSGD: imports GEN4 RespirAct data >
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