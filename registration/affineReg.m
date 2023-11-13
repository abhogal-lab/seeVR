% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <affineReg: uses elastix to perform affine registration of a functional
% to an anatomical image. The function returns the path to both the forward
% and inverse transformation parameter file >
% 
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
% moveImg: path (i.e. directory AND filename with extention) to moving, or
% functional image
%
% moveMask: path (i.e. directory AND filename with extention) to moving, or
% functional mask (binary image)
%
% refImg: path (i.e. directory AND filename with extention) to fixed, or
% anatomical image
%
% refMask: path (i.e. directory AND filename with extention) to fixed, or
% anatomical mask (binary image)
%
% param_af: path to the affine transformation parameter file stored in:
% /seeVR-main/registration/elastix/parameter_files
% this text file can be optimized as needed depending on the data. See the
% elastix manual for details on how this can be done (or contact me)
%
% regdir: path to outputs of affine registration should be saved
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.resultsdir and opts.elastixdir
%
% opts.elastixdir: !!!Important: This is where the elastix binaries are
% stored: e.g., /seeVR-main/registration/elastix

function [image, forward_transform, inverse_transform] = affineReg(moveImg, moveMask, refImg, refMask, param_af, regdir, elastixDir)

if ispc
    elastixDir = fullfile(elastixDir,'windows');
else if ismac
        elastixDir = fullfile(elastixDir,'mac','bin');
    else
        elastixDir = fullfile(elastixDir,'linux','bin');
    end
end

affine_command = [fullfile(elastixDir,'elastix'),' -f ',refImg,' -fmask ',refMask, ' -m ',moveImg,' -mmask ',moveMask,' -p ',param_af,' -out ',regdir ];
dos(affine_command);


forward_transform = fullfile(regdir,'TransformParameters.0.txt')
disp(['transformation parameter file saved as: ',forward_transform])

% calculate inverse 
invdir = fullfile(regdir, 'inverse'); mkdir(invdir);
moveImg_new = fullfile(regdir, 'result.0.nii.gz'); 
inverse_command = [fullfile(elastixDir,'elastix'),' -f ',moveImg,' -fMask ',moveMask, ' -m ',moveImg_new,' -mMask ',refMask,' -t0 ',forward_transform,' -p ',param_af,' -out ',invdir ];
system(inverse_command);

inverse_transform = fullfile(invdir,'TransformParameters.0.txt')
disp(['inverse transformation parameter file saved as: ',inverse_transform])

% load registered image
[image,~] = loadImage(regdir, 'result.0.nii.gz'); 

end

