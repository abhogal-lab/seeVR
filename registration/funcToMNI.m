% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <funcToMNI: uses elastix to register functional image to MNI image
% via the associated structural image >
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

function [] = funcToMNI(moveImg_BET, moveMask, refImg_BET, refImg, refMask, opts)
global opts
tic
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
    elastixrootOS = fullfile(elastixroot,'mac','bin');
else
    elastixrootOS = fullfile(elastixroot,'linux','bin');
end

disp(['moving image: ',moveImg_BET])
if exist(moveImg_BET) == 2
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

disp(['path to reference image: ',refImg_BET])
if exist(refImg_BET) == 2
    disp('brain extracted fixed/reference image found')
else
    error('check brain extracted fixed/reference image filename and extension')
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

%perform first step registration of functional to anatomical image

% affine transformdir
opts.affine_dir = fullfile(opts.resultsdir,'funcToStruct'); mkdir(opts.affine_dir);
disp('performing registration of input to reference image')

%check affine parameter file
disp('checking for affine parameter file for funcToStruct...')

if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileAf.txt')) == 2
    disp('found parameter file')
    param_af_base0 = fullfile(opts.elastixdir,'parameter_files','ParameterFileAf.txt');
    param_af = fullfile(opts.affine_dir,'ParameterFileAf.txt');
    copyfile(param_af_base0, param_af)
    disp(['copying affine parameter file to: ', opts.affine_dir])
else
    error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileAf.txt')])
end

tic
%perform initial func to struct using brain extracted images
[~, forward_transform_step1, inverse_transform_step1] = affineReg(moveImg_BET, moveMask, refImg_BET, refMask, param_af, opts.affine_dir, elastixrootOS);

opts.affineTxParamFile = forward_transform_step1;
opts.inverseAffineTxParamFile = inverse_transform_step1;

name1 = fullfile(opts.affine_dir,'result.0.nii.gz');
name2 = fullfile(opts.affine_dir,'inputToTarget.nii.gz');
movefile(name1, name2)

disp('...')
disp('Registration of functional image to structural image space complete')
disp('Registering structural image to MNI space')

%for non-linear registration, use non-brain extracted images
[trans_params] = structToMNI(refImg, refMask, opts);
disp('registration time....')
toc
opts.regFuncOutputMNI = fullfile(opts.resultsdir,'funcToMNI'); mkdir(opts.regFuncOutputMNI);

opts.affineTxParamFile2 = trans_params.affine_to_target;
opts.bsplineTxParamFile = trans_params.bspline_to_target;
opts.inverseBsplineTxParamFile = trans_params.bspline_to_input;
opts.inverseAffineTxParamFile2 = trans_params.affine_to_input;

mtrans1 = fullfile(opts.affine_dir, 'mTransformParameters.0.txt');
mtrans2 = fullfile(opts.resultsdir,'structToMNI','mTransformParameters.0.txt');
mtrans3 = fullfile(opts.resultsdir,'structToMNI','mTransformParameters.1.txt');

mInvtrans1 = fullfile(opts.regFuncOutputMNI,'mInverseBSTransformParameters.1.txt');
mInvtrans2 = fullfile(opts.regFuncOutputMNI,'mInverseAFTransformParameters.0.txt');
mInvtrans3 = fullfile(opts.regFuncOutputMNI,'mInverseAFTransformParameters.1.txt');

if ispc
    adaptElastixTransFile( opts.affineTxParamFile, mtrans1, 'FinalBSplineInterpolationOrder', '1');
    adaptElastixTransFile( opts.affineTxParamFile2, mtrans2, 'InitialTransformParametersFileName', mtrans1);
    adaptElastixTransFile( mtrans2, mtrans2, 'FinalBSplineInterpolationOrder', '0');
    adaptElastixTransFile( opts.bsplineTxParamFile, mtrans3, 'InitialTransformParametersFileName', mtrans2);
    adaptElastixTransFile( mtrans3, mtrans3, 'FinalBSplineInterpolationOrder', '0');
    %inverse
    adaptElastixTransFile( opts.inverseBsplineTxParamFile, mInvtrans1, 'FinalBSplineInterpolationOrder', '0');
    adaptElastixTransFile( opts.inverseBsplineTxParamFile, mInvtrans1, 'InitialTransformParametersFileName', 'NoInitialTransform');
    adaptElastixTransFile( opts.inverseAffineTxParamFile2, mInvtrans2, 'FinalBSplineInterpolationOrder', '0');

    adaptElastixTransFile( opts.inverseAffineTxParamFile, mInvtrans3, 'InitialTransformParametersFileName', mInvtrans2);
    adaptElastixTransFile( mInvtrans2, mInvtrans2, 'InitialTransformParametersFileName', mInvtrans1);
    
else
    adaptElastixTransFile_linux( opts.affineTxParamFile, mtrans1, 'FinalBSplineInterpolationOrder', '1');
    adaptElastixTransFile_linux( opts.affineTxParamFile2, mtrans2, 'InitialTransformParametersFileName', mtrans1);
    adaptElastixTransFile_linux( mtrans2, mtrans2, 'FinalBSplineInterpolationOrder', '0');
    adaptElastixTransFile_linux( opts.bsplineTxParamFile, mtrans3, 'InitialTransformParametersFileName', mtrans2);
    adaptElastixTransFile_linux( mtrans3, mInvtrans1, 'FinalBSplineInterpolationOrder', '0');

    %inverse
    adaptElastixTransFile_linux( opts.inverseBsplineTxParamFile, mInvtrans1, 'FinalBSplineInterpolationOrder', '0');
    adaptElastixTransFile_linux( opts.inverseAffineTxParamFile2, mInvtrans2, 'FinalBSplineInterpolationOrder', '0');
    adaptElastixTransFile_linux( opts.inverseBsplineTxParamFile, mInvtrans1, 'InitialTransformParametersFileName', 'NoInitialTransform');

    adaptElastixTransFile_linux( opts.inverseAffineTxParamFile, mInvtrans3, 'InitialTransformParametersFileName', mInvtrans2);
    adaptElastixTransFile_linux( mInvtrans2, mInvtrans2, 'InitialTransformParametersFileName', mInvtrans1);
end

%% transform functional image to MNI space
%filename generated by elastix
fileToRename = fullfile(opts.regFuncOutputMNI,'result.nii.gz');

if ispc
        trans_command = [fullfile(elastixrootOS,'transformix'),' -in ',moveImg_BET,' -out ',opts.regFuncOutputMNI,' -tp ',mtrans3 ];
    else
        trans_command = ['transformix -in ',moveImg_BET,' -out ',opts.regFuncOutputMNI,' -tp ',mtrans3];
end
system(trans_command);
 name1 = fileToRename;
 name2 = fullfile(opts.regFuncOutputMNI,'inputToTarget.nii.gz');
 movefile(name1, name2);

%% inverse transform MNI images to functional space
% to mask maps
if ispc
    adaptElastixTransFile( mInvtrans3, mInvtrans3, 'FinalBSplineInterpolationOrder', '0')
else
    adaptElastixTransFile_linux( mInvtrans3, mInvtrans3, 'FinalBSplineInterpolationOrder', '0')
end

maskdir = fullfile(opts.elastixdir,'MNI','labels'); cd(maskdir);
maskname = dir('*.nii.gz*');

fileToRename = fullfile(opts.regFuncOutputMNI,'result.nii.gz');

for kk=1:size(maskname,1)
    maskImg = maskname(kk).name
    [FILEPATH,NAME,EXT] = fileparts(maskImg);
    if ispc
        trans_command = [fullfile(elastixrootOS,'transformix'),' -in ',maskImg,' -out ',opts.regFuncOutputMNI,' -tp ',mInvtrans3 ];
    else
        trans_command = ['transformix -in ',maskImg,' -out ',opts.regFuncOutputMNI,' -tp ',mInvtrans3 ];
    end
    system(trans_command);
    name1 = fileToRename;
    name2 = fullfile(opts.regFuncOutputMNI,[NAME(1:1:end-4),'_toFunc.nii.gz']);
    movefile(name1, name2);
    %move labels files
    name1 = fullfile(maskdir,[maskImg(1:1:end-7),'.txt']);
    name2 = fullfile(opts.regFuncOutputMNI,[maskImg(1:1:end-7),'_labels.txt']);
    if exist(name1) == 2
        copyfile(name1, name2);
    end
end

clear maskname
% to probability maps
if ispc
    adaptElastixTransFile( mInvtrans3, mInvtrans3, 'FinalBSplineInterpolationOrder', '2')
else
    adaptElastixTransFile_linux( mInvtrans3, mInvtrans3, 'FinalBSplineInterpolationOrder', '2')
end

probmaskdir = fullfile(opts.elastixdir,'MNI','prob'); cd(probmaskdir);
maskname = dir('*.nii.gz*');

for kk=1:size(maskname,1)
    maskImg = maskname(kk).name
    [FILEPATH,NAME,EXT] = fileparts(maskImg);
    if ispc
        trans_command = [fullfile(elastixrootOS,'transformix'),' -in ',maskImg,' -out ',opts.regFuncOutputMNI,' -tp ',mInvtrans3 ];
    else
        trans_command = ['transformix -in ',maskImg,' -out ',opts.regFuncOutputMNI,' -tp ',mInvtrans3 ];
    end
    system(trans_command);
    name1 = fileToRename;
    name2 = fullfile(opts.regFuncOutputMNI,[NAME(1:1:end-4),'_toFunc.nii.gz']);
    movefile(name1, name2)
end

toc
end
