function [] = transformMapsStruct(mapdir, TxFile, opts)
global opts

elastixroot = opts.elastixdir;

% Setup OS-dependent paths
if ispc
    elastixrootOS = fullfile(elastixroot, 'windows');
elseif ismac
    elastixrootOS = fullfile(elastixroot, 'mac', 'bin');
else
    elastixrootOS = fullfile(elastixroot, 'linux', 'bin');
end

disp(['Path to parameter maps: ', mapdir])
dirinfo = dir(mapdir);

opts.regOutput = fullfile(mapdir, 'regOutput');
mkdir(opts.regOutput);
disp(['Transformed parameter maps will be saved in: ', opts.regOutput])

% Apply transformix only to NIfTI files
for kk = 1:length(dirinfo)
    if dirinfo(kk).isdir
        continue;
    end

    [~, ~, ext] = fileparts(dirinfo(kk).name);

    % Check for .nii or .nii.gz extension
    isNiiGz = endsWith(dirinfo(kk).name, '.nii.gz');
    isNii = strcmp(ext, '.nii');

    if isNii || isNiiGz
        TxImg = fullfile(dirinfo(kk).folder, dirinfo(kk).name);
        try
            [~] = transformixReg(TxImg, TxFile, opts.regOutput, elastixrootOS);
        catch ME
            warning('Failed to run transformixReg on %s: %s', TxImg, ME.message);
        end
    else
        fprintf('Skipping unsupported file: %s\n', dirinfo(kk).name);
    end
end
