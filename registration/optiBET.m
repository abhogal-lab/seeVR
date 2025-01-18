% Copyright (C) Alex A. Bhogal, 2023, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <optiBET: uses elastix to perform affine + bspline registration of anatomical
% to MNI image. The inverse is then applied and the MNI brain mask is used
% to extract anatomical brain image.
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
% anatImg: full path (i.e. directory AND filename with extention) to anatomical image
%
% opts using T2w image, set opts.T1 = 0;
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.resultsdir and opts.elastixdir
%
% opts.elastixdir: !!!Important: This is where the elastix binaries are
% stored: e.g., /seeVR/registration/elastix

function [anatMask anatBET] = optiBET(anatImg, opts)
global opts;

[anatdir,~,~] = fileparts(anatImg);

try opts.elastixdir; catch
    error('elastix directory not specified for affineReg function... specify path to elastix: opts.elastixdir = ADDPATH')
end

elastixroot = opts.elastixdir;

if isfield(opts,'T1'); else; opts.T1 = 1; end

%setup OS-dependent paths
if ispc
    elastixrootOS = fullfile(elastixroot,'windows');
elseif ismac
    elastixrootOS = fullfile(elastixroot,'mac','bin');
else
    elastixrootOS = fullfile(elastixroot,'linux','bin');
end

if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_BET.txt')) == 2
    disp('found parameter file')
    param_af_base = fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_BET.txt');
    param_af = fullfile(anatdir,'ParameterFileAf.txt');
    copyfile(param_af_base, param_af)
    disp(['copying affine parameter file for forward transform to: ', anatdir])
else
    error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_BET.txt')])
end

if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileBs_BET.txt')) == 2
    disp('found parameter file')
    param_bs_base = fullfile(opts.elastixdir,'parameter_files','ParameterFileBs_BET.txt');
    param_bs = fullfile(anatdir,'ParameterFileBs.txt');
    copyfile(param_bs_base, param_bs)
    disp(['copying affine parameter file for forward transform to: ', anatdir])
else
    error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileBs.txt')])
end

%select MNI image

if opts.T1
    disp('performing initial affine registration of brain extracted moving to brain extracted T1-weighted MNI image')
    refImg = fullfile(opts.elastixdir,'MNI','T1_1mm.nii.gz')
else
    disp('performing initial affine registration of brain extracted moving to brain extracted T2-weighted MNI image')
    refImg = fullfile(opts.elastixdir,'MNI','T2_1mm.nii.gz')
end


if ispc
    bspline_command = [fullfile(elastixrootOS,'elastix'),' -f ',refImg,' -m ',anatImg,' -p ',param_af,' -p ',param_bs,' -out ',anatdir ];

else
    bspline_command = ['elastix -f ',refImg,' -m0 ',anatImg,' -p ',param_af,' -p ',param_bs,' -out ',anatdir ];

end

dos(bspline_command);

opts.affineTxParamFileToTarget = fullfile(anatdir,'TransformParameters.0.txt')
opts.bsplineTxParamFileToTarget = fullfile(anatdir,'TransformParameters.1.txt')

disp(['affinetransformation file saved as: ',fullfile(anatdir,'TransformParameters.0.txt')])
disp(['bsplinetransformation file saved as: ',fullfile(anatdir,'TransformParameters.1.txt')])

%rename output files
if exist(fullfile(anatdir,'result.1.nii.gz')) == 2
    name1 = fullfile(anatdir,'result.1.nii.gz');
    name2 = fullfile(anatdir,'StructToMNI.nii.gz');
    movefile(name1, name2)
end

trans_params = struct();
trans_params.affine_to_target = opts.affineTxParamFileToTarget;
trans_params.bspline_to_target = opts.bsplineTxParamFileToTarget;


outputdir = fullfile(anatdir,'Inverse');
if exist(outputdir) == 7
    cd(outputdir);
    delete *.*
else
    mkdir(outputdir);
end

input_img = fullfile(anatdir,'StructToMNI.nii.gz');

if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_rev_BET.txt')) == 2
    disp('found parameter file')
    param_af_base = fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_rev_BET.txt');
    param_af_rev = fullfile(outputdir,'ParameterFileAf_rev_BET.txt');
    copyfile(param_af_base, param_af_rev)
    disp(['copying affine parameter file for reverse transform to: ', outputdir])
else
    error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_rev.txt')])
end

if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileBs_rev.txt')) == 2
    disp('found parameter file')
    param_bs_base = fullfile(opts.elastixdir,'parameter_files','ParameterFileBs_rev_BET.txt');
    param_bs_rev = fullfile(outputdir,'ParameterFileBs_rev_BET.txt');
    copyfile(param_bs_base, param_bs_rev)
    disp(['copying affine parameter file for reverse transform to: ', outputdir])
else
    error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileBs_rev.txt')])
end

if ispc
    reverse_command = [fullfile(elastixrootOS,'elastix'),' -f ',anatImg,' -m ',input_img,' -p ',param_af_rev,' -p ',param_bs_rev,' -out ',outputdir];
else
    reverse_command = ['elastix -f ',anatImg,' -m ',input_img,' -p ',param_af_rev,' -p ',param_bs_rev,' -out ',outputdir];
end
dos(reverse_command);

opts.affineTxParamFileToInput = fullfile(outputdir,'TransformParameters.0.txt')
opts.bsplineTxParamFileToInput = fullfile(outputdir,'TransformParameters.1.txt')

%rename output files
if exist(fullfile(outputdir,'result.1.nii.gz')) == 2
    name1 = fullfile(outputdir,'result.1.nii.gz');
    name2 = fullfile(outputdir,'invertedT1.nii.gz');
    movefile(name1, name2)
end

trans_params.affine_to_input = opts.affineTxParamFileToInput;
trans_params.bspline_to_input = opts.bsplineTxParamFileToInput;

%apply inverse transform to MNI brain mask
disp('...inverting MNI brain mask')
maskImg= fullfile(opts.elastixdir,'MNI','brain_mask_1mm.nii.gz')

mbs = fullfile(outputdir,'mTransformParameters.1.txt');

if ispc
    adaptElastixTransFile( opts.bsplineTxParamFileToInput, mbs, 'InitialTransformParametersFileName', opts.affineTxParamFileToInput)
    adaptElastixTransFile( mbs, mbs, 'FinalBSplineInterpolationOrder', '0') % nearest neighbor interpolation for binary masks
    trans_command = [fullfile(elastixrootOS,'transformix'),' -in ',maskImg,' -out ',anatdir,' -tp ',mbs ];
else
    adaptElastixTransFile_linux( opts.bsplineTxParamFileToInput, mbs, 'InitialTransformParametersFileName', opts.affineTxParamFileToInput)
    adaptElastixTransFile_linux( mbs, mbs, 'FinalBSplineInterpolationOrder', '0') % nearest neighbor interpolation for binary masks
    trans_command = ['transformix -in ',maskImg,' -out ',anatdir,' -tp ',mbs ];
end

system(trans_command);

%load anatomy and apply mask

[mask,info] = loadMask(anatdir, 'result.nii.gz');
[FILEPATH,NAME,EXT] = fileparts(anatImg);
[anat,anat_info] = loadImage(FILEPATH, [NAME,EXT]);
anat_brain = anat.*(cast(mask, class(anat)));

%save images
niftiwrite(anat_brain,fullfile(anatdir,'anat_brain'),anat_info, 'compressed',1);
anatBET = fullfile(anatdir,'anat_brain.nii.gz')

name1 = fullfile(anatdir, 'result.nii.gz');
name2 = fullfile(anatdir,'anat_brain_mask.nii.gz');
movefile(name1, name2);
anatMask = fullfile(anatdir,'anat_brain_mask.nii.gz')

%clean up files
txtFiles = dir(fullfile(anatdir, '*.txt'));

% Loop through each file and delete it
for i = 1:length(txtFiles)
    filePath = fullfile(anatdir, txtFiles(i).name);
    delete(filePath);
    fprintf('Deleted: %s\n', txtFiles(i).name); % Optional: Print the deleted file name
end

if exist(fullfile(anatdir,'Inverse')) == 7
    status =  rmdir(fullfile(anatdir,'Inverse'), 's');
end

disp('clean-up files completed');
disp('finished brain extraction...');

end

