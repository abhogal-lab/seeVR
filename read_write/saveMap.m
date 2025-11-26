%% Helper function for saving
function saveMap(data, savedir, name, info, opts)

% Handle missing or empty 'info' → fall back to opts.info.map
if nargin < 4 || isempty(info)
    if isfield(opts, 'info') && isfield(opts.info, 'map') && ~isempty(opts.info.map)
        info = opts.info.map;
    else
        error('NIfTI writing requires a valid "info" structure or opts.info.map.');
    end
end

% We now always save via niftiwrite; loadDataSeeVR sets opts.niiwrite = 1
% (keep the field for backwards compatibility but ignore its value)
% opts.niiwrite is assumed == 1 here.

% -----------------------------------------------------------------
% Enforce unit scaling (should already be 1 from loadDataSeeVR, but safe)
% -----------------------------------------------------------------
if isfield(info, 'MultiplicativeScaling')
    info.MultiplicativeScaling = 1;
end
if isfield(info, 'AdditiveOffset')
    info.AdditiveOffset = 0;
end
if isfield(info, 'raw')
    if isfield(info.raw, 'scl_slope')
        info.raw.scl_slope = 1;
    end
    if isfield(info.raw, 'scl_inter')
        info.raw.scl_inter = 0;
    end
end

% -----------------------------------------------------------------
% Make sure info.Datatype matches data class and cast data accordingly
% -----------------------------------------------------------------
info = syncInfoDatatypeWithData(info, class(data));
data = castDataToDatatype(data, info.Datatype);

% -----------------------------------------------------------------
% Write compressed NIfTI
% -----------------------------------------------------------------
niftiwrite(data, fullfile(savedir, name), info, "Compressed", true);

end

% ================== LOCAL HELPERS ================== %

function info = syncInfoDatatypeWithData(info, dataClass)
    switch dataClass
        case 'double'
            dtStr = 'double';
        case 'single'
            dtStr = 'single';
        case 'int16'
            dtStr = 'int16';
        case 'uint16'
            dtStr = 'uint16';
        case 'uint8'
            dtStr = 'uint8';
        case 'int8'
            dtStr = 'int8';
        case 'int32'
            dtStr = 'int32';
        case 'uint32'
            dtStr = 'uint32';
        otherwise
            % Fallback: if info already has a datatype, keep it; else use double
            if isfield(info, 'Datatype') && ~isempty(info.Datatype)
                dtStr = info.Datatype;
            else
                dtStr = 'double';
            end
    end
    info.Datatype = dtStr;
end

function dataOut = castDataToDatatype(dataIn, dtStr)
    switch dtStr
        case 'double'
            dataOut = double(dataIn);
        case 'single'
            dataOut = single(dataIn);
        case 'int16'
            dataOut = int16(dataIn);
        case 'uint16'
            dataOut = uint16(dataIn);
        case 'uint8'
            dataOut = uint8(dataIn);
        case 'int8'
            dataOut = int8(dataIn);
        case 'int32'
            dataOut = int32(dataIn);
        case 'uint32'
            dataOut = uint32(dataIn);
        otherwise
            % Unknown / unsupported string → leave as is
            dataOut = dataIn;
    end
end
