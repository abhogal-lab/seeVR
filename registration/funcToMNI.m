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

function [] = funcToMNI(moveImg, moveMask, refImg, refMask, mapdir, opts)
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

%perform first step registration of functional to anatomical image

% affine transformdir
opts.affine_dir = fullfile(opts.resultsdir,'funcToStruct'); mkdir(opts.affine_dir);
disp('performing registration of input to reference image')

[~, forward_transform_step1, inverse_transform_step1] = affineReg(moveImg, moveMask, refImg, refMask, param_af, opts.affine_dir, elastixrootOS);

opts.affineTxParamFile = forward_transform_step1;
opts.inverseAffineTxParamFile = inverse_transform_step1;

name1 = fullfile(opts.affine_dir,'result.0.nii.gz');
name2 = fullfile(opts.affine_dir,'funcToStruct.nii.gz');
movefile(name1, name2)

%apply transforms
disp(['path to parameter maps: ',mapdir])
dirinfo = dir(mapdir);

opts.regOutput = fullfile(mapdir,'funcToStruct'); mkdir(opts.regOutput);
disp(['transformed parameter maps will be saved in: ', opts.regOutput])

%apply Tx
for kk=1:size(dirinfo,1)
    if dirinfo(kk).isdir
    else
        TxImg = fullfile(dirinfo(kk).folder,dirinfo(kk).name);
        [~] = transformixReg(TxImg, opts.affineTxParamFile, opts.regOutput, elastixrootOS);
    end
end

disp('...')
disp('Registration of functional image to structural image space complete')
disp('Associated parameter maps transformed to structural image space')
disp('Registering structural image to MNI space')

[trans_params] = structToMNI(refImg, refMask, opts);

opts.regOutputMNI = fullfile(mapdir,'funcToMNI'); mkdir(opts.regOutputMNI);

trans1 = opts.affineTxParamFile;
trans2 = trans_params.affine_to_target;
trans3 = trans_params.bspline_to_target;
mtrans1 = fullfile(opts.affine_dir, 'mTransformParameters.0.txt');
mtrans2 = fullfile(opts.resultsdir,'structToMNI','mTransformParameters.0.txt');
mtrans3 = fullfile(opts.resultsdir,'structToMNI','mTransformParameters.1.txt');

if ispc
    adaptElastixTransFile( trans1, mtrans1, 'FinalBSplineInterpolationOrder', '0');
    adaptElastixTransFile( trans2, mtrans2, 'InitialTransformParametersFileName', mtrans1);
    adaptElastixTransFile( mtrans2, mtrans2, 'FinalBSplineInterpolationOrder', '0');
    adaptElastixTransFile( trans3, mtrans3, 'InitialTransformParametersFileName', mtrans2);
    adaptElastixTransFile( mtrans3, mtrans3, 'FinalBSplineInterpolationOrder', '1');
else
    adaptElastixTransFile_linux( trans1, mtrans1, 'FinalBSplineInterpolationOrder', '0');
    adaptElastixTransFile_linux( trans2, mtrans2, 'InitialTransformParametersFileName', mtrans1);
    adaptElastixTransFile_linux( trans2, mtrans2, 'FinalBSplineInterpolationOrder', '0');
    adaptElastixTransFile_linux( trans3, mtrans3, 'InitialTransformParametersFileName', mtrans2);
    adaptElastixTransFile_linux( mtrans3, mtrans3, 'FinalBSplineInterpolationOrder', '1');
end

%apply Tx
for kk=1:size(dirinfo,1)
    if dirinfo(kk).isdir
    else
        TxImg = fullfile(dirinfo(kk).folder,dirinfo(kk).name);
        [~] = transformixReg(TxImg, mtrans3, opts.regOutputMNI, elastixrootOS);
    end
end

disp('...')
disp('Registration of functional image to MNI space complete')
disp('Associated parameter maps transformed to MNI space')

toc
end
