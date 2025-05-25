%% Helper function for saving
function saveMap(data, savedir, name, info, opts)

% Handle missing or empty 'info'
if nargin < 4 || isempty(info)
    info = [];
end

% Write NIfTI or save as image
if opts.niiwrite
    if isempty(info)
        error('NIfTI writing requires a valid "info" structure.');
    end
    niftiwrite(data, fullfile(savedir, name), info, "Compressed", 1);
else
    saveImageData(data, opts.headers.map, savedir, name, 64);
end

end
