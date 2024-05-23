% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <funcToStruct: uses elastix to register a functional image to an anatomical image >
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
% opts: options structure containing required variables for this specific
% function; i.e. opts.resultsdir and opts.elastixdir
%
% opts.elastixdir: !!!Important: This is where the elastix binaries are
% stored: e.g., /seeVR-main/registration/elastix
%
% Parameter files are stored in
% /seeVR-main/registration/elastix/parameter_files and can be optimized as
% needed. See the elastix manual

function [] = funcToStruct(moveImg, moveMask, refImg, refMask, opts)
global opts
try opts.elastixdir; catch
    error('elastix directory not specified... specify OS-dependent path to elastix: opts.elastixdir = ADDPATH')
end

disp(['moving image: ',moveImg])
if exist(moveImg) == 2
    disp('moving image found')
else
    error('check moving image filename and extension')
end

disp(['moving mask: ',moveMask])
if exist(moveMask) == 2
    disp('moving mask found')
else
    error('check moving mask filename and extension')
end

disp(['path to reference image: ',refImg])
if exist(refImg) == 2
    disp('fixed/reference image found')
else
    error('check fixed/reference image filename and extension')
end

disp(['path to reference image: ',refMask])
if exist(refMask) == 2
    disp('fixed/reference mask found')
else
    error('check fixed/reference mask filename and extension')
end

%perform affine registration of func to anat

% affine transformdir
opts.affine_dir = fullfile(opts.resultsdir,'affineReg'); mkdir(opts.affine_dir);
disp(['moving to fixed transform saved in: ',opts.affine_dir]);

disp('checking for parameter file...')

if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileAf.txt')) == 2
    disp('found parameter file')
    param_af = fullfile(opts.elastixdir,'parameter_files','ParameterFileAf.txt');
else
    error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileAf.txt')])
end

disp('performing registration of moving to fixed image')

[~, forward_transform, ~] = affineReg(moveImg, moveMask, refImg, refMask, param_af, opts.affine_dir, opts.elastixdir);

name1 = fullfile(opts.affine_dir,'result.0.nii.gz');
name2 = fullfile(opts.affine_dir,'func2anat.nii.gz');
movefile(name1, name2)

disp('registration of moving to anat complete')
disp(['transformation file saved as: ',fullfile(opts.affine_dir,'TransformParameters.0.txt')])
disp(['final image saved as: ',name2])
opts.affineTxParamFile = fullfile(opts.affine_dir,'TransformParameters.0.txt');

end
