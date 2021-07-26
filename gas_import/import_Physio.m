%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

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