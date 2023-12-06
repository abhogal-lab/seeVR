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
global opts

if ispc
    elastixDir = fullfile(opts.elastixDir,'windows');
else if ismac
        elastixDir = fullfile(opts.elastixDir,'mac','bin');
    else
        elastixDir = fullfile(opts.elastixDir,'linux','bin');
    end
end

transform_command = [fullfile(opts.elastixDir, 'transformix'),' -in ',inputImg,' -out ',transformdir, ' -tp ', transformfile];
dos(transform_command);

% load registered image
[image,~] = loadImage(transformdir, 'result.nii.gz'); 

end