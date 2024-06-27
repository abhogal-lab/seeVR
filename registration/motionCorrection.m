function [data_mcf] = motionCorrection(data,opts)

global opts
try opts.elastixdir; catch
    error('elastix directory not specified... specify OS-dependent path to elastix: opts.elastixdir = ADDPATH')
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

if exist(fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_mcf.txt')) == 2
    disp('found parameter file')
    param_af = fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_mcf.txt');
else
    error(['check elastix parameter file. Expected: ',fullfile(opts.elastixdir,'parameter_files','ParameterFileAf_mcf.txt')])
end
opts.outputdir_mc = fullfile(opts.savedir,'motion_corrected');
mkdir(opts.outputdir_mc);

reference = squeeze(data(:,:,:,1));

cd(opts.outputdir_mc)
niftiwrite(cast(reference, opts.mapDatatype),'reference_image',opts.info.map);
ref_img = fullfile(opts.outputdir_mc,'reference_image.nii');

data_mcf = zeros(size(data));
data_tr = zeros(size(data));

affine_transmat = fullfile(opts.outputdir_mc,'TransformParameters.0.txt');
outputdir_tr = fullfile(opts.outputdir_mc,'transformation_files');
mkdir(outputdir_tr);

for ii=1:size(data,4)

    disp(['registering volume number ',int2str(ii)]);

    cd(opts.outputdir_mc)
    moving_img = squeeze(data(:,:,:,ii));
    niftiwrite(cast(moving_img, opts.mapDatatype),'moving_image',opts.info.map);

    moving = fullfile(opts.outputdir_mc,'moving_image.nii');

    if ispc
        affine_command = [fullfile(elastixrootOS,'elastix'),' -f ',ref_img,' -m ',moving,' -p ',param_af,' -out ',opts.outputdir_mc,' > ouput.txt'];
    else
        affine_command = ['elastix -f ',ref_img,' -m ',moving,' -p ',param_af,' -out ',opts.outputdir_mc,' > ouput.txt'];
    end

    dos(affine_command);

    %load registered image
    [img, ~] = loadImage(opts.outputdir_mc,'result.0.nii.gz');
    data_mcf(:,:,:,ii) = img;

    %%
    if ispc
    rename_command = ['copy ',affine_transmat,' ',outputdir_tr,'\TransformParameters.',int2str(ii),'.txt'];
    else
    rename_command = ['cp ',affine_transmat,' ',outputdir_tr,'\TransformParameters.',int2str(ii),'.txt'];
    end

    dos(rename_command);

end
disp('...completed motion correction')
cd(opts.outputdir_mc)
data_mcf = cast(data_mcf, opts.tsDatatype);
niftiwrite(data_mcf,'data_mcf',opts.info.ts);

end
