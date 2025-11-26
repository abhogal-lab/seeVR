function [image, info, header] = loadDataSeeVR(pathname, filename, varargin)
%LOADDATASEEVR Unified loader for NIfTI timeseries, masks, maps, anatomical images.
%
%   [image, info, header] = loadDataSeeVR(pathname, filename, 'type', TYPE, 'useUntouch', FLAG)
%
%   TYPE:
%     'auto'       - default, treat as anatomical/image
%     'timeseries' - fMRI / 4D data, initializes opts.func* and TR/dyn
%     'ts'         - alias for 'timeseries'
%     'mask'       - binary/label mask
%     'map'        - parametric map
%
%   useUntouch:
%     [] (default) - automatically uses load_untouch_nii for int16 data
%     true/false   - force use or non-use of load_untouch_nii
%
%   Design:
%   - Always uses niftiinfo to obtain MATLAB-style INFO structs.
%   - Optionally uses load_untouch_nii just for voxel data when needed.
%   - Forces all scaling (MultiplicativeScaling / scl_slope / scl_inter) to
%     1 (and additive offsets to 0) but stores original values in
%     opts.originalScaling.* so they are not lost.
%   - For timeseries: if TR cannot be determined, prompts user via inputdlg
%     to enter TR in seconds.
%
%   This makes opts.info.* always consistent with MATLAB niftiwrite/niftiinfo
%   expectations, so saving can be unified through niftiwrite.

warning('off');
global opts;
if isempty(opts)
    opts = struct;
end

% -------------------- INPUT PARSING -------------------- %
p = inputParser;
p.addParameter('type', 'auto');
p.addParameter('useUntouch', []);
p.parse(varargin{:});
loadType   = lower(p.Results.type);
useUntouch = p.Results.useUntouch;

% -------------------- BASIC PATH & INFO ---------------- %
% Some systems append newline at end of filename
if ~ispc && ~isempty(filename) && filename(end) == char(10)
    filename(end) = [];
end

image_path = fullfile(pathname, filename);

% Always get MATLAB-style info
info = niftiinfo(image_path);

% Decide if we want untouch for voxel data (e.g. int16 quirks)
if isempty(useUntouch)
    useUntouch = strcmpi(info.Datatype, 'int16');
end

% -------------------- LOAD VOXEL DATA ------------------ %
header = [];
if useUntouch
    % Use untouch only for voxel array + raw header
    nii    = load_untouch_nii(image_path);
    image  = double(nii.img);
    header = nii.hdr;
else
    % Standard MATLAB NIfTI reader
    image  = niftiread(image_path);
    header = [];  % not needed for saving, but kept for backwards compatibility
end

% If info somehow failed (very rare), derive minimal structure
if isempty(info) && ~isempty(header) && isfield(header, 'dime')
    info = struct( ...
        'PixelDimensions', header.dime.pixdim(2:5), ...
        'ImageSize',       header.dime.dim(2:5),    ...
        'Datatype',        header.dime.datatype);
end

% Clean infinities
image(isinf(image)) = 0;

% Voxel dimensions (first 3 elements)
voxDim = info.PixelDimensions(1:min(3, numel(info.PixelDimensions)));

% ==============================================================
%   CAPTURE ORIGINAL SCALING + FORCE UNIT SCALING
% ==============================================================

origScale = struct( ...
    'info_MultiplicativeScaling', [], ...
    'info_AdditiveOffset',        [], ...
    'raw_scl_slope',              [], ...
    'raw_scl_inter',              [], ...
    'header_scl_slope',           [], ...
    'header_scl_inter',           []);

% MATLAB info-layer scaling
if isfield(info, 'MultiplicativeScaling')
    origScale.info_MultiplicativeScaling = info.MultiplicativeScaling;
    info.MultiplicativeScaling = 1;
end
if isfield(info, 'AdditiveOffset')
    origScale.info_AdditiveOffset = info.AdditiveOffset;
    info.AdditiveOffset = 0;
end

% Raw NIfTI scaling (scl_slope / scl_inter) inside info.raw
if isfield(info, 'raw')
    if isfield(info.raw, 'scl_slope')
        origScale.raw_scl_slope = info.raw.scl_slope;
        info.raw.scl_slope      = 1;
    end
    if isfield(info.raw, 'scl_inter')
        origScale.raw_scl_inter = info.raw.scl_inter;
        info.raw.scl_inter      = 0;
    end
end

% Header scaling (only meaningful if we used untouch)
if ~isempty(header) && isfield(header, 'dime')
    if isfield(header.dime, 'scl_slope')
        origScale.header_scl_slope = header.dime.scl_slope;
        header.dime.scl_slope      = 1;
    end
    if isfield(header.dime, 'scl_inter')
        origScale.header_scl_inter = header.dime.scl_inter;
        header.dime.scl_inter      = 0;
    end
end

% Ensure opts.originalScaling exists
if ~isfield(opts, 'originalScaling') || isempty(opts.originalScaling)
    opts.originalScaling = struct;
end

% NOTE: we now *always* intend to save via niftiwrite → no need to branch
opts.niiwrite = 1;

% ==============================================================
%   TYPE-SPECIFIC PROCESSING → ALWAYS USING MATLAB-STYLE INFO
% ==============================================================

switch loadType
    %% ----------------- TIMESERIES ----------------- %%
    case {'timeseries','ts'}
        opts.funcpath  = pathname;
        opts.funcfile  = filename;
        opts.voxelsize = voxDim;

        hasTR = false;
        trVal = [];

        if ndims(image) > 3
            % Try to read TR (pixdim(5)) from info.raw or untouch header
            if isfield(info, 'raw') && isfield(info.raw, 'pixdim') && numel(info.raw.pixdim) >= 5 ...
                    && ~isempty(info.raw.pixdim(5))
                trVal = info.raw.pixdim(5);
            elseif ~isempty(header) && isfield(header, 'dime') && numel(header.dime.pixdim) >= 5
                trVal = header.dime.pixdim(5);
            end

            if ~isempty(trVal)
                % Convert ms → s if it looks like ms
                if trVal > 10
                    trVal = trVal / 1000;
                end
                opts.TR = trVal;
                hasTR   = true;
            end

            opts.dyn = size(image, 4);
        end

        % If TR not found → request from user
        if ~hasTR
            answer = inputdlg({'Enter TR (seconds):'}, 'TR Required', 1, {''});
            if isempty(answer)
                error('Timeseries detected but TR not provided — operation cancelled by user.');
            end
            userTR = str2double(answer{1});
            if isnan(userTR) || userTR <= 0
                error('Invalid TR value. TR must be a positive numeric value in seconds.');
            end

            opts.TR = userTR;
            fprintf('User-specified TR = %.3f seconds\n', userTR);
        end

        % Store info in MATLAB-native format
        opts.info.ts    = info;
        opts.tsDatatype = info.Datatype;

        % Initialize map info from timeseries (same structure as niftiinfo)
        opts.info.map = info;
        if isfield(opts.info.map, 'raw')
            opts.info.map.raw.dim(1)    = 3;
            opts.info.map.raw.dim(5)    = 1;
            opts.info.map.raw.pixdim(5) = 0;
        end
        if isfield(opts.info.map, 'PixelDimensions')
            opts.info.map.PixelDimensions = opts.info.ts.PixelDimensions(1:3);
        end
        if isfield(opts.info.map, 'ImageSize')
            opts.info.map.ImageSize = opts.info.ts.ImageSize(1:3);
        end
        opts.mapDatatype = info.Datatype;

        % Store original scaling for ts/map/mask derived from this
        opts.originalScaling.ts   = origScale;
        opts.originalScaling.map  = origScale;
        opts.originalScaling.mask = origScale;

    %% ----------------- MASK ----------------- %%
    case 'mask'
        opts.info.mask = info;
        if isfield(opts.info.mask, 'MultiplicativeScaling')
            opts.info.mask.MultiplicativeScaling = 1;
        end
        opts.maskDatatype = info.Datatype;

        opts.originalScaling.mask = origScale;

    %% ----------------- MAP ----------------- %%
    case 'map'
        opts.info.map    = info;
        opts.mapDatatype = info.Datatype;

        opts.originalScaling.map = origScale;

    %% ----------------- IMAGE / ANATOMY / AUTO ------- %%
    otherwise  % {'auto','image','anat'}
        opts.voxelsize_image = voxDim;
        opts.imagepath       = pathname;
        opts.imagefile       = filename;
        opts.info.image      = info;
        opts.imageDatatype   = info.Datatype;

        opts.info.anat    = info;
        opts.anatDatatype = info.Datatype;

        opts.originalScaling.image = origScale;
        opts.originalScaling.anat  = origScale;
end

end
