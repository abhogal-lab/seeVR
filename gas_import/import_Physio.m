% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <import_Physio: imports GEN4 RespirAct data >
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

function [MRTimes2, ID, FRCmL, VdmL, TissuestoreO2mL, TissuestoreCO2mL, VO2mLmin, VCO2mLmin, QmLmin, hBconcentrationgdLBlood, Restingmetabolicscalefactor, ResponseReason] = import_Physio(filename, dataLines)


%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 12);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["MRTimes2", "ID", "FRCmL", "VdmL", "TissuestoreO2mL", "TissuestoreCO2mL", "VO2mLmin", "VCO2mLmin", "QmLmin", "hBconcentrationgdLBlood", "Restingmetabolicscalefactor", "ResponseReason"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "ResponseReason", "TrimNonNumeric", true);
opts = setvaropts(opts, "ResponseReason", "ThousandsSeparator", ",");

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
MRTimes2 = tbl.MRTimes2;
ID = tbl.ID;
FRCmL = tbl.FRCmL;
VdmL = tbl.VdmL;
TissuestoreO2mL = tbl.TissuestoreO2mL;
TissuestoreCO2mL = tbl.TissuestoreCO2mL;
VO2mLmin = tbl.VO2mLmin;
VCO2mLmin = tbl.VCO2mLmin;
QmLmin = tbl.QmLmin;
hBconcentrationgdLBlood = tbl.hBconcentrationgdLBlood;
Restingmetabolicscalefactor = tbl.Restingmetabolicscalefactor;
ResponseReason = tbl.ResponseReason;
end