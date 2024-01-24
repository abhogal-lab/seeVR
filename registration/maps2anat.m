% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <maps2anat: uses elastix to register paramter maps an anatomical image 
% via the associated functional image >
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
% mapdir: path to where parameter files (derived from functional data) are
% stored. These should have the same dimensions as the functional data
% (except for -t dimention). All .nii* files in this directory will be
% transformed to anat space
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.resultsdir and opts.elastixdir
%
% opts.elastixdir: !!!Important: This is where the elastix binaries are
% stored: e.g., /seeVR-main/registration/elastix
%
% parameter files are stored in
% /seeVR-main/registration/elastix/parameter_files and can be optimized as
% needed. See the elastix manual

function [] = maps2anat(moveImg, moveMask, refImg, refMask, mapdir, opts)
global opts

try opts.elastixdir; catch
    error('elastix directory not specified... specify path to elastix: opts.elastixdir = ADDPATH')
end

elastixroot = opts.elastixdir;
elastixparam = fullfile(elastixroot, 'parameter_files')

disp('checking for parameter file...')

if exist(fullfile(elastixparam,'ParameterFileAf.txt')) == 2
    disp('found parameter file')
    param_af = fullfile(elastixparam,'ParameterFileAf.txt');
else 
    error(['check elastix parameter file. Expected: ',fullfile(elastixparam,'ParameterFileAf.txt')])
end

%setup OS-dependent paths
if ispc
elastixrootOS = fullfile(elastixroot,'windows');
elseif ismac
elastixrootOS = fullfile(elastixroot,'mac');
else
elastixrootOS = fullfile(elastixroot,'linux');
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

disp('performing registration of moving to fixed image')

[~, forward_transform, inverse_transform] = affineReg(moveImg, moveMask, refImg, refMask, param_af, opts.affine_dir, elastixrootOS);
opts.affineTxParamFile = forward_transform;
opts.inverseAffineTxParamFile = inverse_transform;

name1 = fullfile(opts.affine_dir,'result.0.nii.gz');
name2 = [fullfile(opts.affine_dir,'func2anat.nii.gz')];
movefile(name1, name2)

%apply transforms
disp(['path to parameter maps: ',mapdir])
dirinfo = dir(mapdir);

opts.regOutput = fullfile(mapdir,'regOutput'); mkdir(opts.regOutput);
disp(['transformed parameter maps will be saved in: ', opts.regOutput])

%apply Tx
for kk=1:size(dirinfo,1)
    if dirinfo(kk).isdir
    else
        TxImg = fullfile(dirinfo(kk).folder,dirinfo(kk).name);
        [~] = transformixReg(TxImg, opts.affineTxParamFile, opts.regOutput, elastixrootOS);
        %rename result image
        [FILEPATH,NAME,EXT] = fileparts(TxImg);
        name1 = fullfile(opts.regOutput,'result.nii.gz');
        name2 = fullfile(opts.regOutput,[NAME,'_anat.nii.gz']);
        movefile(name1, name2)
    end
end

disp('...')
disp('registration of func to anat complete')
disp('parameter maps transformed to anat')
end
