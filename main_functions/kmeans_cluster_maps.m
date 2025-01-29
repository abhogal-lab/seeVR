% Copyright (C) Alex A. Bhogal, 2025, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <cluster_maps: performs k-means clusting using input maps
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% *************************************************************************
%
%% maps: Cell array containing 3D parameter maps {map1, map2, ..., mapN}
%
% numClusters: Number of clusters to segment the volume into
%
% clusteredVolume: 3D volume with voxels labeled by cluster index
%
% clusterCenters: Centers of the clusters


function [clusteredVolume, clusterCenters] = cluster_maps(maps, numClusters, opts)

% Get dimensions of the 3D maps
global opts

[xDim, yDim, zDim] = size(maps{1});
numMaps = length(maps);

% Reshape the 3D maps into a 2D matrix (each row is a voxel, each column is a feature)
featureMatrix = zeros(xDim * yDim * zDim, numMaps);
for i = 1:numMaps
    featureMatrix(:, i) = reshape(maps{i}, [], 1);
end

% Remove any NaN or Inf values that may corrupt clustering
validIdx = all(~isnan(featureMatrix), 2) & all(~isinf(featureMatrix), 2);
validFeatures = featureMatrix(validIdx, :);

% Perform K-means clustering on the valid features
[clusterIdx, clusterCenters] = kmeans(validFeatures, numClusters, 'Replicates', 10, 'MaxIter', 200);

% Initialize the output volume with NaNs
clusteredVolume = NaN(xDim * yDim * zDim, 1);
clusteredVolume(validIdx) = clusterIdx;

% Reshape clusteredVolume back to 3D
clusteredVolume = reshape(clusteredVolume, [xDim, yDim, zDim]);

if opts.niiwrite
    cd(opts.resultsdir);
    niftiwrite(cast(clusteredVolume, opts.mapDatatype),['kmClusterMap_',int2str(numClusters),'_clusters'],opts.info.map);
else
    saveImageData(clusteredVolume, opts.headers.map, opts.resultsdir, ['kmClusterMap_',int2str(numClusters),'_clusters.nii.gz'], 64);
end
end