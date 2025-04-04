%% Helper function for saving
function saveMap(data, name, info, opts)
if opts.niiwrite
    niftiwrite(data, fullfile(opts.ALFFdir, name), info);
else
    saveImageData(data, opts.headers.map, opts.ALFFdir, name, 64);
end
end