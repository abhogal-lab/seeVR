function [clusteredVolume, clusterCenters] = cluster_maps(maps, numClusters, opts)
global opts
% Copyright (C) Alex A. Bhogal, 2025, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <cluster_maps: performs k-means clustering using input maps

% -------------------------------------------------------------------------
% opts.gpu = 1  -> attempt GPU usage
% opts.gpu = 0  -> CPU only
%
% Notes:
% - GPU execution for kmeans depends on MATLAB/toolbox/version support.
% - If GPU kmeans is not supported, the function falls back to CPU.
% -------------------------------------------------------------------------

if nargin < 3 || isempty(opts)
    opts = struct;
end
if ~isfield(opts, 'gpu') || isempty(opts.gpu)
    opts.gpu = 0;
end

useGPU = logical(opts.gpu);

% Get dimensions of the 3D maps
[xDim, yDim, zDim] = size(maps{1});
numMaps = length(maps);

% Check that all maps have identical dimensions
for i = 2:numMaps
    if ~isequal(size(maps{i}), [xDim, yDim, zDim])
        error('All maps must have identical dimensions.');
    end
end

% Choose numeric class based on first map
outClass = class(maps{1});
if ~isfloat(maps{1})
    outClass = 'single';
end

% ---------------------------------------------------------------------
% Build feature matrix
% ---------------------------------------------------------------------
if useGPU
    try
        g = gpuDevice; %#ok<NASGU>
        featureMatrix = gpuArray.zeros(xDim * yDim * zDim, numMaps, outClass);

        for i = 1:numMaps
            featureMatrix(:, i) = reshape(gpuArray(cast(maps{i}, outClass)), [], 1);
        end

        validIdx = all(~isnan(featureMatrix), 2) & all(~isinf(featureMatrix), 2);
        validFeatures = featureMatrix(validIdx, :);

        % Try GPU kmeans
        try
            [clusterIdx, clusterCenters] = kmeans(validFeatures, numClusters, ...
                'Replicates', 10, 'MaxIter', 200);

            % Gather outputs back to CPU
            clusterIdx = gather(clusterIdx);
            clusterCenters = gather(clusterCenters);
            validIdx = gather(validIdx);

        catch
            warning('GPU kmeans failed or is unsupported in this MATLAB version. Falling back to CPU');

            % Fall back to CPU
            featureMatrix = gather(featureMatrix);
            validIdx = gather(validIdx);
            validFeatures = featureMatrix(validIdx, :);

            [clusterIdx, clusterCenters] = kmeans(validFeatures, numClusters, ...
                'Replicates', 10, 'MaxIter', 200);
        end

    catch
        warning('GPU requested but unavailable. Falling back to CPU');

        % CPU path
        featureMatrix = zeros(xDim * yDim * zDim, numMaps, outClass);
        for i = 1:numMaps
            featureMatrix(:, i) = reshape(cast(maps{i}, outClass), [], 1);
        end

        validIdx = all(~isnan(featureMatrix), 2) & all(~isinf(featureMatrix), 2);
        validFeatures = featureMatrix(validIdx, :);

        [clusterIdx, clusterCenters] = kmeans(validFeatures, numClusters, ...
            'Replicates', 10, 'MaxIter', 200);
    end

else
    % -----------------------------------------------------------------
    % CPU path
    % -----------------------------------------------------------------
    featureMatrix = zeros(xDim * yDim * zDim, numMaps, outClass);
    for i = 1:numMaps
        featureMatrix(:, i) = reshape(cast(maps{i}, outClass), [], 1);
    end

    validIdx = all(~isnan(featureMatrix), 2) & all(~isinf(featureMatrix), 2);
    validFeatures = featureMatrix(validIdx, :);

    [clusterIdx, clusterCenters] = kmeans(validFeatures, numClusters, ...
        'Replicates', 10, 'MaxIter', 200);
end

% ---------------------------------------------------------------------
% Rebuild clustered volume
% ---------------------------------------------------------------------
clusteredVolume = NaN(xDim * yDim * zDim, 1);
clusteredVolume(validIdx) = clusterIdx;
clusteredVolume = reshape(clusteredVolume, [xDim, yDim, zDim]);

% Save output
% Original header
infostruct = opts.info.map;

% New size from clustered output
newSize = size(clusteredVolume);

% Old size
oldSize = infostruct.ImageSize;

% Update ImageSize
infostruct.ImageSize = newSize;

% --- Adjust voxel size to preserve FOV ---
if isfield(infostruct, 'PixelDimensions') && numel(oldSize) >= 3
    scaleFactor = oldSize(1:3) ./ newSize(1:3);
    infostruct.PixelDimensions(1:3) = infostruct.PixelDimensions(1:3) .* scaleFactor;
end

% Optional: ensure datatype matches
infostruct.Datatype = class(clusteredVolume);
infostruct.BitsPerPixel = 8 * whos('clusteredVolume').bytes / numel(clusteredVolume);

% Save
saveMap(clusteredVolume, opts.resultsdir, ...
    ['gmmClusterMap_', int2str(numClusters), '_clusters'], ...
    infostruct, opts);

end