% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <transformixReg: is used to apply elastix transformation parameter files
% to, for example, register parameter maps to anatomical or atlas space
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
% *************************************************************************
% Registration relies on elastix. For more information see the Elastix
% hompage at: https://elastix.lumc.nl/
%
% inputImg: path (i.e. directory AND filename with extention) to moving, or
% functional image
%
% moveMask: path (i.e. directory AND filename with extention) to moving
% image that is to be transformed - e.g., parameter map
%
% transformfile: transformation file to be applied. This can be the output
% of the affineReg.m function
%
% transformdir: directory to save transformed images
%
% opts.elastixdir: !!!Important: This is where the elastix binaries are
% stored: e.g., /seeVR-main/registration/elastix

function [image] = transformixReg(inputImg, transformfile, transformdir, opts)
global opts;

try opts.elastixdir; catch
    error('elastix directory not specified fir affineReg function... specify path to elastix: opts.elastixdir = ADDPATH')
end

elastixroot = opts.elastixdir;

%setup OS-dependent paths
if ispc
    elastixrootOS = fullfile(elastixroot,'windows');
elseif ismac
    elastixrootOS = fullfile(elastixroot,'mac','bin');
else
    elastixrootOS = fullfile(elastixroot,'linux','bin');
end

if ispc
    transform_command = [fullfile(elastixrootOS, 'transformix'),' -in ',inputImg,' -out ',transformdir, ' -tp ', transformfile];
else
    transform_command = ['transformix -in ',inputImg,' -out ',transformdir, ' -tp ', transformfile];
end
dos(transform_command);

%rename image

[~, name, ext] = fileparts(inputImg);

name1 = fullfile(transformdir,'result.nii.gz');
name2 = fullfile(transformdir,['tr_',name,'.nii.gz']);
movefile(name1, name2);
% load registered image
[image,~] = loadImage(transformdir, ['tr_',name,'.nii.gz']);

end