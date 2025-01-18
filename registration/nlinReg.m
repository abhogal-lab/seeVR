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
% stored: e.g., /seeVR/registration/elastix

function [trans_params] = nlinReg(moveImg, moveMask, refImg, refMask, regdir, opts)
global opts;

try opts.elastixdir; catch
    error('elastix directory not specified for affineReg function... specify path to elastix: opts.elastixdir = ADDPATH')
end

if isfield(opts,'invert_bspline'); else; opts.invert_bspline = 1; end

elastixroot = opts.elastixdir;

%setup OS-dependent paths
if ispc
    elastixrootOS = fullfile(elastixroot,'windows');
elseif ismac
    elastixrootOS = fullfile(elastixroot,'mac','bin');
else
    elastixrootOS = fullfile(elastixroot,'linux','bin');
end

if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileAf.txt')) == 2
    disp('found parameter file')
    param_af_base = fullfile(opts.elastixdir,'parameter_files','ParameterFileAf.txt');
    param_af = fullfile(regdir,'ParameterFileAf_rev.txt');
    copyfile(param_af_base, param_af)
    disp(['copying affine parameter file for forward transform to: ', regdir])
else
    error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileAf.txt')])
end

if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileBs.txt')) == 2
    disp('found parameter file')
    param_bs_base = fullfile(opts.elastixdir,'parameter_files','ParameterFileBs_rev.txt');
    param_bs = fullfile(regdir,'ParameterFileBs_rev.txt');
    copyfile(param_bs_base, param_bs)
    disp(['copying affine parameter file for forward transform to: ', regdir])
else
    error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileBs.txt')])
end

if ispc
    bspline_command = [fullfile(elastixrootOS,'elastix'),' -f ',refImg,' -m ',moveImg,' -p ',param_af,' -p ',param_bs,' -out ',regdir ];

else
    bspline_command = ['elastix -f ',refImg,' -m0 ',moveImg,' -p ',param_af,' -p ',param_bs,' -out ',regdir ];

end

dos(bspline_command);

opts.affineTxParamFileToTarget = fullfile(regdir,'TransformParameters.0.txt')
opts.bsplineTxParamFileToTarget = fullfile(regdir,'TransformParameters.1.txt')

disp(['affinetransformation file saved as: ',fullfile(regdir,'TransformParameters.0.txt')])
disp(['bsplinetransformation file saved as: ',fullfile(regdir,'TransformParameters.1.txt')])

%rename output files
if exist(fullfile(opts.bspline_dir,'result.1.nii.gz')) == 2
    name1 = fullfile(opts.bspline_dir,'result.1.nii.gz');
    name2 = fullfile(opts.bspline_dir,'InputToTarget_nLin.nii.gz');
    movefile(name1, name2)
end

trans_params = struct();

trans_params.affine_to_target = opts.affineTxParamFileToTarget;
trans_params.bspline_to_target = opts.bsplineTxParamFileToTarget;

if opts.invert_bspline

    outputdir = fullfile(regdir,'Inverse');
    if exist(outputdir) == 7
        cd(outputdir);
        delete *.*
    else
        mkdir(outputdir);
    end

    input_img = fullfile(opts.bspline_dir,'InputToTarget_nLin.nii.gz');

    if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_rev.txt')) == 2
        disp('found parameter file')
        param_af_base = fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_rev.txt');
        param_af_rev = fullfile(outputdir,'ParameterFileAf_rev.txt');
        copyfile(param_af_base, param_af_rev)
        disp(['copying affine parameter file for reverse transform to: ', outputdir])
    else
        error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_rev.txt')])
    end

    if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileBs_rev.txt')) == 2
        disp('found parameter file')
        param_bs_base = fullfile(opts.elastixdir,'parameter_files','ParameterFileBs_rev.txt');
        param_bs_rev = fullfile(outputdir,'ParameterFileBs_rev.txt');
        copyfile(param_bs_base, param_bs_rev)
        disp(['copying affine parameter file for reverse transform to: ', outputdir])
    else
        error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileBs_rev.txt')])
    end

    if ispc
        reverse_command = [fullfile(elastixrootOS,'elastix'),' -f ',moveImg,' -m ',input_img,' -p ',param_af_rev,' -p ',param_bs_rev,' -out ',outputdir];
    else
        reverse_command = ['elastix -f ',moveImg,' -m ',input_img,' -p ',param_af_rev,' -p ',param_bs_rev,' -out ',outputdir];
    end
    dos(reverse_command);

    opts.affineTxParamFileToInput = fullfile(outputdir,'TransformParameters.0.txt')
    opts.bsplineTxParamFileToInput = fullfile(outputdir,'TransformParameters.1.txt')

    %rename output files
    if exist(fullfile(outputdir,'result.1.nii.gz')) == 2
        name1 = fullfile(outputdir,'result.1.nii.gz');
        name2 = fullfile(outputdir,'TargetToInput_nLin.nii.gz');
        movefile(name1, name2)
    end

    trans_params.affine_to_input = opts.affineTxParamFileToInput;
    trans_params.bspline_to_input = opts.bsplineTxParamFileToInput;
end



end

