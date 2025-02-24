function [] = plotMap(sourceImg, mask, paramMap, map, filename, foldername, rotations, opts)
% plotMap: Plots parameter map next to source image in transverse, coronal, and sagittal views.
%          Slices are rotated first, then cropped to a 2D bounding box. The figure window is sized
%          based on the number of rows/columns to keep the subplots close together.
%
% USAGE:
%   plotMap(sourceImg, mask, paramMap, map, filename, foldername, rotations, opts)
%
% INPUTS:
%   sourceImg   : [X x Y x Z (x T)] source image (e.g., mean BOLD or anatomical).
%   mask        : [X x Y x Z] binary mask defining voxels of interest.
%   paramMap    : [X x Y x Z (x T)] parameter map of same size as sourceImg (3D or 4D).
%   map         : colormap for the paramMap slices.
%   filename    : base filename for saving images (no extension).
%   foldername  : path to directory where images will be saved.
%   rotations   : 1x3 numeric vector specifying # of 90-degree rotations (2D) for each orientation:
%                   rotations(1) -> transverse
%                   rotations(2) -> coronal
%                   rotations(3) -> sagittal
%                 e.g., [0 1 -1] => no rotation for transverse,
%                                    rotate coronal 90° CCW,
%                                    rotate sagittal 90° clockwise.
%   opts        : structure with optional fields:
%       opts.scale  - data range (default: [-5 5])
%       opts.row    - number of subplot rows (default: 5)
%       opts.col    - number of subplot columns (default: 6)
%       opts.step   - step size through slices (default: 1)
%       opts.start  - starting slice index (orientation-dependent default)
%
% OUTPUT:
%   Saves PNG and EPS images for each orientation in 'foldername', named:
%       filename_ori1.png / filename_ori1.eps   (transverse)
%       filename_ori2.png / filename_ori2.eps   (coronal)
%       filename_ori3.png / filename_ori3.eps   (sagittal)
%
% -------------------------------------------------------------------------
global opts;  % If your workflow relies on global opts

% ------------------- Set Default opts Values -------------------
if ~isfield(opts, 'scale'), opts.scale = [-5 5]; end
if ~isfield(opts, 'row'),   opts.row   = 5;     end
if ~isfield(opts, 'col'),   opts.col   = 6;     end
if ~isfield(opts, 'step'),  opts.step  = 1;     end

% ------------------- Validate foldername -------------------
if ~exist('foldername','var') || isempty(foldername) || ~ischar(foldername)
    foldername = pwd;  % fallback
end

% ------------------- Validate rotations -------------------
if ~exist('rotations','var') || isempty(rotations) || numel(rotations) < 3
    rotations = [0 0 0];  % default = no rotation
end

% ------------------- Convert Inputs to Double -------------------
sourceImg = double(sourceImg);
mask      = double(mask); 
mask(mask > 0) = 1;  % ensure strictly 0/1
paramMap  = double(paramMap);

% Dimensions
dims = size(mask);
% orientation: 1=transverse (vary 3rd dim)
%              2=coronal    (vary 2nd dim)
%              3=sagittal   (vary 1st dim)
defaultStart = [floor(dims(3)/3), floor(dims(2)/3), floor(dims(1)/3)];

% ------------------- Main Loop Over 3 Orientations -------------------
for ori = 1:3
    % Determine the start slice for this orientation
    if ~isfield(opts,'start') || isempty(opts.start)
        startSlice = defaultStart(ori);
    else
        startSlice = opts.start;
    end
    % Clamp startSlice within dimension range
    switch ori
        case 1, startSlice = min(max(startSlice,1), dims(3)); % transverse
        case 2, startSlice = min(max(startSlice,1), dims(2)); % coronal
        case 3, startSlice = min(max(startSlice,1), dims(1)); % sagittal
    end
    
    % Prepare to store slices for bounding box determination
    nPairs = floor((opts.row * opts.col)/2);  % each pair => (source, paramMap)
    sourceSlices = cell(nPairs,1);
    paramSlices  = cell(nPairs,1);
    actualCount  = 0;  % how many valid slices we gather
    
    % Gather slices for bounding box
    for p = 1:nPairs
        jj       = 2*p - 1;  % 1,3,5,... for pairs
        sliceIdx = startSlice + jj * opts.step;
        
        % Check dimension bounds
        switch ori
            case 1, if sliceIdx<1 || sliceIdx>dims(3), continue; end
            case 2, if sliceIdx<1 || sliceIdx>dims(2), continue; end
            case 3, if sliceIdx<1 || sliceIdx>dims(1), continue; end
        end
        
        sSlice = getSlice(sourceImg, mask, sliceIdx, ori);
        pSlice = getSlice(paramMap,  mask, sliceIdx, ori);
        pSlice = nanmean(pSlice,4);  % if 4D
        
        % Outside mask => 0 => convert to NaN for alpha
        sSlice(sSlice == 0) = NaN;
        pSlice(pSlice == 0) = NaN;
        
        % Rotate first
        kRot = rotations(ori);
        if kRot ~= 0
            sSlice = rot90(sSlice, kRot);
            pSlice = rot90(pSlice, kRot);
        end
        
        actualCount            = actualCount + 1;
        sourceSlices{actualCount} = sSlice;
        paramSlices{actualCount}  = pSlice;
    end
    
    if actualCount < 1
        warning('No valid slices for orientation %d. Skipping...', ori);
        continue;
    end
    
    % Find global bounding box of all slices (both source & param)
    [rMin, rMax, cMin, cMax] = findGlobalBoundingBox(sourceSlices(1:actualCount), ...
                                                     paramSlices(1:actualCount));
    
    % Create a figure sized by row/col to keep them close together.
    % Adjust this scaleFactor to shrink or enlarge your subplots.
    scaleFactor = 140;  % e.g. 120 pixels per tile
    figWidth  = opts.col * scaleFactor;
    figHeight = opts.row * scaleFactor;
    
    f = figure('Units','pixels',...
               'Position',[100, 100, figWidth, figHeight],...
               'Color','k');
    
    t = tiledlayout(opts.row, opts.col, 'TileSpacing','none', 'Padding','none');
    
    pairIndex = 1;
    for p = 1:nPairs
        if pairIndex > actualCount, break; end
        
        sSlice = sourceSlices{pairIndex};
        pSlice = paramSlices{pairIndex};
        
        % Crop to the bounding box
        sCropped = sSlice(rMin:rMax, cMin:cMax);
        pCropped = pSlice(rMin:rMax, cMin:cMax);
        
        % ---------- Source Tile ----------
        nexttile(t);
        h1 = imagesc(sCropped);
        colormap(gray); freezeColors;
        axis image off;
        set(h1,'AlphaData', ~isnan(sCropped));  % transparent outside mask
        
        % ---------- ParamMap Tile ----------
        nexttile(t);
        h2 = imagesc(pCropped, opts.scale);
        colormap(map); freezeColors;
        axis image off;
        set(h2,'AlphaData', ~isnan(pCropped));
        
        pairIndex = pairIndex + 1;
    end
    
    drawnow;  % finalize layout
    
    % Save figure
    outBase = fullfile(foldername, [filename, '_ori', num2str(ori)]);
    exportgraphics(f,[outBase, '.png'],...
        'Resolution',600, 'BackgroundColor','black');
    exportgraphics(f,[outBase, '.eps'],...
        'Resolution',600, 'BackgroundColor','black');
    
    close(f);
end
end % plotMap

% ------------------- Helper: getSlice -------------------
function slice2D = getSlice(vol3D, mask3D, sliceIndex, orientation)
% Extract a 2D slice (plus any extra 4th-dim if present) according to:
%   ori=1 (transverse): vol3D(:,:,sliceIndex)
%   ori=2 (coronal)   : vol3D(:,sliceIndex,:)
%   ori=3 (sagittal)  : vol3D(sliceIndex,:,:)
% Then multiply by mask3D so that outside the mask is zero. We'll convert
% zeros to NaNs for alpha transparency in the main code.

switch orientation
    case 1 % Transverse
        if ndims(vol3D) == 3
            slice2D = mask3D(:,:,sliceIndex) .* vol3D(:,:,sliceIndex);
        else
            slice2D = mask3D(:,:,sliceIndex) .* vol3D(:,:,sliceIndex,:);
        end
        
    case 2 % Coronal
        if ndims(vol3D) == 3
            slice2D = squeeze(mask3D(:,sliceIndex,:) .* vol3D(:,sliceIndex,:));
        else
            slice2D = squeeze(mask3D(:,sliceIndex,:)) .* squeeze(vol3D(:,sliceIndex,:,:));
        end
        
    case 3 % Sagittal
        if ndims(vol3D) == 3
            slice2D = squeeze(mask3D(sliceIndex,:,:) .* vol3D(sliceIndex,:,:));
        else
            slice2D = squeeze(mask3D(sliceIndex,:,:)) .* squeeze(vol3D(sliceIndex,:,:,:));
        end
end

end

% ------------------- Helper: findGlobalBoundingBox -------------------
function [rMin, rMax, cMin, cMax] = findGlobalBoundingBox(sourceSlices, paramSlices)
% Given cell arrays of 2D source slices & param slices, each possibly
% rotated, find the minimal row/column bounding box that includes all non-NaN
% pixels across *all* slices. This ensures consistent cropping.

rMin = Inf; rMax = 0;
cMin = Inf; cMax = 0;

nSlices = numel(sourceSlices);
for i = 1:nSlices
    sS = sourceSlices{i};
    pS = paramSlices{i};
    if ~isequal(size(sS), size(pS))
        error('Source slice and param slice differ in size at index %d.', i);
    end
    
    combined = ~isnan(sS) | ~isnan(pS);
    if ~any(combined(:)), continue; end
    
    [rs, cs] = find(combined);
    rMin = min(rMin, min(rs));
    rMax = max(rMax, max(rs));
    cMin = min(cMin, min(cs));
    cMax = max(cMax, max(cs));
end

% If no valid voxels, default to [1 1 1 1] to avoid errors
if rMin == Inf || rMax == 0
    rMin = 1; rMax = 1;
    cMin = 1; cMax = 1;
end

end
