% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <structToMNI: uses elastix to perform registration of structural image
% to MNI image.
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
% anatomical image. 
%
% moveMask: path (i.e. directory AND filename with extention) to moving, or
% anatomical mask (binary image)
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

function [trans_params] = structToMNI(moveImg, moveMask, opts)
global opts

if isfield(opts,'T1'); else; opts.T1 = 1; end
if isfield(opts,'useBET'); else; opts.useBET = 0; end

try opts.elastixdir; catch
    error('elastix directory not specified... specify OS-dependent path to elastix: e.g. opts.elastixdir = /.../seeVR/registration/elastix/')
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

opts.bspline_dir = fullfile(opts.resultsdir,'structToMNI');
if exist(opts.bspline_dir) == 7
    cd(opts.bspline_dir);
    delete *.*
else
    mkdir(opts.bspline_dir);
end
%select MNI image

if opts.T1
    disp('performing initial affine registration of brain extracted moving to brain extracted T1-weighted MNI image')
    refImg_BET = fullfile(opts.elastixdir,'MNI','T1_1mm_brain.nii.gz')
    refImg = fullfile(opts.elastixdir,'MNI','T1_1mm.nii.gz') 
else
    disp('performing initial affine registration of brain extracted moving to brain extracted T2-weighted MNI image')
    refImg_BET = fullfile(opts.elastixdir,'MNI','T2_1mm_brain.nii.gz')
    refImg = fullfile(opts.elastixdir,'MNI','T2_1mm.nii.gz')
end

refMask = fullfile(opts.elastixdir,'MNI','brain_mask_1mm.nii.gz')

disp('checking for affine parameter file...')

if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileAf.txt')) == 2
    disp('found parameter file')
    param_af_base = fullfile(opts.elastixdir,'parameter_files','ParameterFileAf.txt');
    param_af = fullfile(opts.bspline_dir,'ParameterFileAf.txt');
    copyfile(param_af_base, param_af)
    disp(['copying affine parameter file to: ', opts.bspline_dir])
else
    error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileAf.txt')])
end

disp('checking for bspline parameter file...')

if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileBs.txt')) == 2
    disp('found parameter file')
    param_bs_base = fullfile(opts.elastixdir,'parameter_files','ParameterFileBs.txt');
    param_bs = fullfile(opts.bspline_dir,'ParameterFileBs.txt');
    copyfile(param_bs_base, param_bs)
    disp(['copying affine parameter file to: ', opts.bspline_dir])
else
    error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileBs.txt')])
end
if opts.useBET;  refImg_BET = refImg; end

[trans_params] = nlinReg(moveImg, moveMask, refImg_BET, refMask, opts.bspline_dir, opts.elastixdir);

% update transform parameters
outputdir = fullfile(opts.bspline_dir,'Inverse');
mbs = fullfile(outputdir,'mTransformParameters.1.txt');

if ispc
    adaptElastixTransFile( trans_params.bspline_to_input, mbs, 'InitialTransformParametersFileName', trans_params.affine_to_input)
    adaptElastixTransFile( mbs, mbs, 'FinalBSplineInterpolationOrder', '0') % nearest neighbor interpolation for binary masks
else
    adaptElastixTransFile_linux( trans_params.bspline_to_input, mbs, 'InitialTransformParametersFileName', trans_params.affine_to_input)
    adaptElastixTransFile_linux( mbs, mbs, 'FinalBSplineInterpolationOrder', '0') % nearest neighbor interpolation for binary masks
end

%% apply transformations to MNI masks

maskdir = fullfile(opts.elastixdir,'MNI','labels'); cd(maskdir);
maskname = dir('*.nii.gz*');

fileToRename = fullfile(outputdir,'result.nii.gz');

for kk=1:size(maskname,1)
    maskImg = maskname(kk).name
    [FILEPATH,NAME,EXT] = fileparts(maskImg);
    if ispc
        trans_command = [fullfile(elastixrootOS,'transformix'),' -in ',maskImg,' -out ',outputdir,' -tp ',mbs ];
    else
        trans_command = ['transformix -in ',maskImg,' -out ',outputdir,' -tp ',mbs ];
    end
    system(trans_command);
    name1 = fileToRename;
    name2 = fullfile(outputdir,[NAME(1:1:end-4),'_toInput.nii.gz']);
    movefile(name1, name2);
    %move labels files
    name1 = fullfile(maskdir,[maskImg(1:1:end-7),'.txt']);
    name2 = fullfile(outputdir,[maskImg(1:1:end-7),'_labels.txt']);
    if exist(name1) == 2
        copyfile(name1, name2);
    end
end

clear maskname
% to probability maps
if ispc
    adaptElastixTransFile( mbs, mbs, 'FinalBSplineInterpolationOrder', '2') 
else
    adaptElastixTransFile_linux( mbs, mbs, 'FinalBSplineInterpolationOrder', '2') 
end

probmaskdir = fullfile(opts.elastixdir,'MNI','prob'); cd(probmaskdir);
maskname = dir('*.nii.gz*');

for kk=1:size(maskname,1)
    maskImg = maskname(kk).name
    [FILEPATH,NAME,EXT] = fileparts(maskImg);
    if ispc
        trans_command = [fullfile(elastixrootOS,'transformix'),' -in ',maskImg,' -out ',outputdir,' -tp ',mbs ];
    else
        trans_command = ['transformix -in ',maskImg,' -out ',outputdir,' -tp ',mbs ];
    end
    system(trans_command);
    name1 = fileToRename;
    name2 = fullfile(outputdir,[NAME(1:1:end-4),'_toInput.nii.gz']);
    movefile(name1, name2)
end
end

