function [map, voxelCount, meanVal] = atlasMaps(data, atlas)
% generate atlas maps

atlasMap =  zeros(size(atlas));
voxelMap =  zeros(size(atlas));
maxInd = max(atlas(:));

for ii=1:maxInd

    mask = atlas;
    mask(mask == ii) = -1;
    mask(mask > 0) = 0; mask = abs(mask);
    tmp = mask.*cast(data, class(mask));
    tmp = tmp(:); tmp(tmp == 0) = NaN;

    voxels = tmp; voxels(isnan(voxels)) = [];
    voxelCount(ii) = numel(voxels);
    meanVal(ii) = nanmean(tmp);

    atlasMap(atlas == ii) = meanVal(ii);
    voxelMap(atlas == ii) = voxelCount(ii);

    disp(['calculating region :', int2str(ii), ' of ',int2str(maxInd)])
end

atlasMap = cast(atlasMap, class(data));
voxelMap = cast(voxelMap, class(data));
map.atlasMap = atlasMap;
map.voxelMap = voxelMap;
voxelCount = voxelCount(:);
meanVal = meanVal(:);
end