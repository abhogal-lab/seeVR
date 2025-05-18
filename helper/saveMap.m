%% Helper function for saving
function saveMap(data, savedir, name, info, opts)
if opts.niiwrite
    niftiwrite(data, fullfile(savedir, name), info, "Compressed",1);
else
    saveImageData(data, opts.headers.map, savedir, name, 64);
end
end