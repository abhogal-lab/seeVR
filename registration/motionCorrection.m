function [data_mcf, motion_params] = motionCorrection(data,opts)

global opts
if ~isfield(opts,'info.map');
    disp('no header present, loading input data to produce headers')
    [pathstr, name, ext] = fileparts(data)
    [data,~] = loadTimeseries(pathstr,[name,ext]);
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

if opts.use_mean
    reference = mean(data,4); %use mean image
else
    reference = squeeze(data(:,:,:,1)); %use first image
end

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
        affine_command = strjoin([fullfile(elastixrootOS,'elastix'),' -f ',ref_img,' -m ',moving,' -p ',param_af,' -out ',opts.outputdir_mc,' > ouput.txt']);
    else
        affine_command = strjoin(['elastix -f ',ref_img,' -m ',moving,' -p ',param_af,' -out ',opts.outputdir_mc,' > ouput.txt']);
    end

    dos(affine_command);

    %load registered image
    [img, ~] = loadImage(opts.outputdir_mc,'result.0.nii.gz');
    data_mcf(:,:,:,ii) = img;

    %%
    if ispc
        rename_command = strjoin(['copy ',affine_transmat,' ',outputdir_tr,'\TransformParameters.',int2str(ii),'.txt']);
    else
        rename_command = strjoin(['cp ',affine_transmat,' ',outputdir_tr,'\TransformParameters.',int2str(ii),'.txt']);
    end

    dos(rename_command);

end

disp('...completed motion correction')
cd(opts.outputdir_mc)
data_mcf = cast(data_mcf, opts.tsDatatype);
niftiwrite(data_mcf,'data_mcf',opts.info.ts);

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
            param_vals = str2num(param_str{1}); 
            motion_params = [motion_params; param_vals];
            break;
        end
        tline = fgetl(fid);
    end
    fclose(fid);
end

output_filename = fullfile(opts.outputdir_mc,'motion_params.txt')
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
