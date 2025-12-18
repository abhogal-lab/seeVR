function [data_mcf, motion_params] = motionCorrection(data,opts)

global opts
if ~isfield(opts,'info.map');
    disp('no header present, loading input data to produce headers')
    [pathstr, name, ext] = fileparts(data);
    [data,~] = loadTimeseries(pathstr,[name,ext]);
else
    [pathstr, ~, ~] = fileparts(opts.info.map.Filename);
end

global opts
try opts.elastixdir; catch
    error('elastix directory not specified... specify OS-dependent path to elastix: opts.elastixdir = ADDPATH')
end

if isfield(opts,'use_mean'); else; opts.use_mean = 1; end

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

% =====================================================================
% Helper: map MATLAB class to NIfTI Datatype and BitsPerPixel
% =====================================================================
inClass = class(data);   % we will preserve this

    function info = syncDatatype(info, cls)
        switch cls
            case 'double'
                info.Datatype     = 'double';
                info.BitsPerPixel = 64;
            case 'single'
                info.Datatype     = 'single';
                info.BitsPerPixel = 32;
            case 'uint8'
                info.Datatype     = 'uint8';
                info.BitsPerPixel = 8;
            case 'int8'
                info.Datatype     = 'int8';
                info.BitsPerPixel = 8;
            case 'uint16'
                info.Datatype     = 'uint16';
                info.BitsPerPixel = 16;
            case 'int16'
                info.Datatype     = 'int16';
                info.BitsPerPixel = 16;
            case 'uint32'
                info.Datatype     = 'uint32';
                info.BitsPerPixel = 32;
            case 'int32'
                info.Datatype     = 'int32';
                info.BitsPerPixel = 32;
            otherwise
                warning('Unsupported input class %s, defaulting to double in NIfTI.', cls);
                info.Datatype     = 'double';
                info.BitsPerPixel = 64;
        end
        % make sure PixelDimensions has enough entries
        if numel(info.PixelDimensions) < numel(info.ImageSize)
            info.PixelDimensions(end+1:numel(info.ImageSize)) = 1;
        end
        % remove low-level raw header if present (avoids headerBytes assert)
        if isfield(info,'raw')
            info = rmfield(info,'raw');
        end
    end

% =====================================================================
% Reference image (only saving logic changed)
% =====================================================================
if opts.use_mean
    reference = mean(data,4); %use mean image
else
    reference = squeeze(data(:,:,:,1)); %use first image
end

cd(opts.outputdir_mc)

refInfo           = opts.info.map;
refInfo.ImageSize = size(reference);
refInfo           = syncDatatype(refInfo, inClass);

niftiwrite(cast(reference, inClass),'reference_image',refInfo,'Compressed',false,'Version','NIfTI1');
ref_img = fullfile(opts.outputdir_mc,'reference_image.nii');

% =====================================================================
% Main loop (Elastix calls unchanged, only moving-image saving fixed)
% =====================================================================
data_mcf = [];
data_tr  = []; %#ok<NASGU>

affine_transmat = fullfile(opts.outputdir_mc,'TransformParameters.0.txt');
outputdir_tr    = fullfile(opts.outputdir_mc,'transformation_files');
mkdir(outputdir_tr);

for ii=1:size(data,4)

    disp(['registering volume number ',int2str(ii)]);

    cd(opts.outputdir_mc)
    moving_img = squeeze(data(:,:,:,ii));

    movInfo           = opts.info.map;
    movInfo.ImageSize = size(moving_img);
    movInfo           = syncDatatype(movInfo, inClass);

    niftiwrite(cast(moving_img, inClass),'moving_image',movInfo,'Compressed',false,'Version','NIfTI1');
    moving = fullfile(opts.outputdir_mc,'moving_image.nii');

    if ispc
        affine_command = ([fullfile(elastixrootOS,'elastix'),' -f ',ref_img,' -m ',moving,' -p ',param_af,' -out ',opts.outputdir_mc,' > ouput.txt']);
    else
        affine_command = (['elastix -f ',ref_img,' -m ',moving,' -p ',param_af,' -out ',opts.outputdir_mc,' > ouput.txt']);
    end

    dos(affine_command);

    %load registered image
    [img, ~] = loadImage(opts.outputdir_mc,'result.0.nii.gz');

    if ii == 1
        % Allocate based on what Elastix actually returns
        inClass = class(img);  % likely 'single' if ResultImagePixelType "float"
        data_mcf = zeros([size(img) size(data,4)], inClass);
    end

    data_mcf(:,:,:,ii) = img;  % no cast needed, types match

    %%
    if ispc
        rename_command = (['copy ',affine_transmat,' ',outputdir_tr,'\TransformParameters.',int2str(ii),'.txt']);
    else
        rename_command = (['cp ',affine_transmat,' ',outputdir_tr,'\TransformParameters.',int2str(ii),'.txt']);
    end

    dos(rename_command);

end

disp('...completed motion correction')

% =====================================================================
% Final save of data_mcf (fixed)
% =====================================================================
cd(opts.outputdir_mc)

if isfield(opts,'info') && isfield(opts.info,'ts')
    outInfo = opts.info.ts;
else
    outInfo = opts.info.map;
end

outInfo.ImageSize = size(data_mcf);
outInfo           = syncDatatype(outInfo, inClass);

niftiwrite(cast(data_mcf, inClass),'data_mcf',outInfo,'Compressed',false,'Version','NIfTI1');

% =====================================================================
% Motion parameters (unchanged)
% =====================================================================
cd(fullfile(opts.outputdir_mc,'transformation_files'))

disp('generating motion file')
names = ls("Transform*");
motion_params = [];

for idx=1:length(names)
    filename = ['TransformParameters.',int2str(idx),'.txt'];
    filename(filename == ' ') = [];
    %extract params

    fid = fopen(filename, 'r');
    % Read lines until finding TransformParameters
    tline = fgetl(fid);
    while ischar(tline)
        if contains(tline, '(TransformParameters')
            % Extract numerical parameters
            param_str = extractBetween(tline, '(TransformParameters ', ')');
            param_vals = str2num(param_str{1}); %#ok<ST2NM>
            motion_params = [motion_params; param_vals];
            break;
        end
        tline = fgetl(fid);
    end
    fclose(fid);
end

output_filename = fullfile(opts.outputdir_mc,'motion_params.txt');
% Write the concatenated parameters to output file
writematrix(motion_params, output_filename, 'Delimiter',  '\t');
disp('motion parameters are saved in motion_params.txt file')

% Plotting each motion parameter
n_params = size(motion_params, 2);
figure('Position', [100, 100, 400, 100*n_params]);

cmap = parula(n_params);
for i = 1:n_params
    subplot(n_params,1,i);
    plot(motion_params(:,i), 'Color', cmap(i,:), 'LineWidth', 1.5);
    ylabel(sprintf('Param #%d', i));
    xlabel('Frame Index');
    grid on;
    xlim([0 length(motion_params)])
end
sgtitle('Motion Parameters Over Time');

saveas(gcf,fullfile(opts.outputdir_mc,'motion_parameters.png'))
saveas(gcf,fullfile(opts.outputdir_mc,'motion_parameters.fig'))
close(gcf)

cd(pathstr);
end
