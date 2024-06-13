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
global opts;

try opts.elastixdir; catch
    error('elastix directory not specified for affineReg function... specify path to elastix: opts.elastixdir = ADDPATH')
end

if isfield(opts,'invert_affine'); else; opts.invert_affine = 1; end

elastixroot = opts.elastixdir;

%setup OS-dependent paths
if ispc
    elastixrootOS = fullfile(elastixroot,'windows');
    affine_command = [fullfile(elastixrootOS,'elastix'),' -f ',refImg,' -fmask ',refMask, ' -m ',moveImg,' -mmask ',moveMask,' -p ',param_af,' -out ',regdir ];
elseif ismac
    elastixrootOS = fullfile(elastixroot,'mac','bin'); %not tested - may be buggy
    affine_command = [fullfile(elastixrootOS,'elastix'),' -f ',refImg,' -fmask ',refMask, ' -m ',moveImg,' -mmask ',moveMask,' -p ',param_af,' -out ',regdir ];
else
    elastixrootOS = fullfile(elastixroot,'linux','bin');
    affine_command = ['elastix -f ',refImg,' -fmask ',refMask, ' -m ',moveImg,' -mmask ',moveMask,' -p ',param_af,' -out ',regdir ];
end

system(affine_command);

forward_transform = fullfile(regdir,'TransformParameters.0.txt')
disp(['transformation parameter file saved as: ',forward_transform])

if opts.invert_affine
    
    % calculate inverse
    invdir = fullfile(regdir, 'inverse'); mkdir(invdir);
    moveImg_new = fullfile(regdir, 'result.0.nii.gz');
    
    if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_rev.txt')) == 2
        disp('found parameter file')
        param_af_base = fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_rev.txt');
        param_af_rev = fullfile(invdir,'ParameterFileAf_rev.txt');
        copyfile(param_af_base, param_af_rev)
        disp(['copying affine parameter file for reverse transform to: ', invdir])
    else
        error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_rev.txt')])
    end
    
    if ispc
        inverse_command = [fullfile(elastixrootOS,'elastix'),' -f ',moveImg, ' -m ',moveImg_new,' -t0 ',forward_transform,' -p ',param_af_rev,' -out ',invdir ];
    else
        inverse_command = ['elastix -f ',moveImg,' -m ',moveImg_new,' -t0 ',forward_transform,' -p ',param_af_rev,' -out ',invdir ];
    end
    
    system(inverse_command);
    
    inverse_transform = fullfile(invdir,'TransformParameters.0.txt')
    disp(['inverse transformation parameter file saved as: ',inverse_transform])
    disp('applying inverse transformation to reference image')
    
    [~] = transformixReg(refImg, inverse_transform, invdir, elastixrootOS);
    
else
    inverse_transform = [];
end

% load registered image
[image,~] = loadImage(regdir, 'result.0.nii.gz');

end

