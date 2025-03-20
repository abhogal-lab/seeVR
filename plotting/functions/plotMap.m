function [] = plotMap(sourceImg, mask, paramMap, map, filename, foldername, rotations, opts)
% plotMap: Plots parameter map next to source image in transverse, coronal, and sagittal views.
%
% USAGE:
%   plotMap(sourceImg, mask, paramMap, map, filename, foldername, rotations, opts)
%
% INPUTS:
%   sourceImg   : source image (e.g., mean BOLD or anatomical image).
%   mask        : binary mask defining voxels of interest.
%   paramMap    : parameter map with the same size as sourceImg (can be 3D or 4D).
%   map         : colormap to be used when displaying the paramMap slices.
%   filename    : base filename for saving images (no extension).
%   foldername  : path to directory where images will be saved.
%   rotations   : 1x3 numeric vector specifying the number of 90-degree rotations
%                 to apply for each orientation:
%                    rotations(1) -> transverse
%                    rotations(2) -> coronal
%                    rotations(3) -> sagittal
%                 e.g., [0 1 -1] means no rotation for transverse, rotate coronal
%                 slices 90° CCW, and rotate sagittal slices 90° clockwise (=-1).
%
%   opts        : structure with optional fields:
%       opts.scale  - data range               (default: [-5 5])
%       opts.row    - number of rows in figure (default: 5)
%       opts.col    - number of columns        (default: 6)
%       opts.step   - step size through slices (default: 1)
%       opts.start  - starting slice index     (orientation-dependent default)
%
% OUTPUT:
%   Saved PNG and EPS images for each orientation in 'foldername', named:
%       filename_ori1.png / filename_ori1.eps   (transverse)
%       filename_ori2.png / filename_ori2.eps   (coronal)
%       filename_ori3.png / filename_ori3.eps   (sagittal)
%
% -------------------------------------------------------------------------
global opts;  % If your workflow relies on global opts

% ================= Default Options =================
if ~isfield(opts, 'scale'), opts.scale = [-5 5]; end
if ~isfield(opts, 'row'),   opts.row   = 5;     end
if ~isfield(opts, 'col'),   opts.col   = 6;     end
if ~isfield(opts, 'step'),  opts.step  = 1;     end
tight = 0;

% ================= Validate foldername =================
if ~exist('foldername','var') || isempty(foldername) || ~ischar(foldername)
    foldername = pwd;  % fallback if invalid
end

% ================= Validate rotations =================
if ~exist('rotations','var') || isempty(rotations) || numel(rotations) < 3
    rotations = [0 0 0];  % default = no rotation
end

% ================= Convert Inputs to Double =================
sourceImg = double(sourceImg);
mask      = double(mask);
mask(mask > 0) = 1;   % ensure mask is strictly 0/1
paramMap  = double(paramMap);

% ================= Dimensions & Default Start Indices =================
dims = size(mask);
% orientation: 1=transverse (vary 3rd dim)
%              2=coronal    (vary 2nd dim)
%              3=sagittal   (vary 1st dim)
defaultStart = [floor(dims(3)/3), floor(dims(2)/3), floor(dims(1)/3)];

% ================= Main Loop over 3 Orientations =================
for ori = 1:3
    
    % Determine start slice for this orientation
    if ~isfield(opts,'start') || isempty(opts.start)
        startSlice = defaultStart(ori);
    else
        startSlice = opts.start;
        % Make sure we don't go out of bounds
        switch ori
            case 1 % transverse
                if startSlice > dims(3), startSlice = floor(dims(3)/2); end
            case 2 % coronal
                if startSlice > dims(2), startSlice = floor(dims(2)/2); end
            case 3 % sagittal
                if startSlice > dims(1), startSlice = floor(dims(1)/2); end
        end
    end
    
    % Create a figure for this orientation
    if tight
    %figure('Units','pixels',...
    %       'Position',[10, 100, 100 * opts.col * size(paramMap,2), ...
    %                   opts.row * size(paramMap,2)],...
    %       'Color','k');
    else
    figure('Color','k');
    end

    % We'll create pairs of subplots: (sourceSlice, paramMapSlice)
    for jj = 1 : 2 : (opts.row*opts.col - 1)
        
        sliceIdx = startSlice + jj * opts.step;

        % ------------------ Source Slice ------------------
        sourceSlice = getSlice(sourceImg, mask, sliceIdx, ori);
        
        % ------------------ ParamMap Slice ------------------
        pmSlice = getSlice(paramMap, mask, sliceIdx, ori);
        % If paramMap is 4D, average across the 4th dimension
        pmSlice = nanmean(pmSlice,4);
        
        % Set zeros to NaN for display
        pmSlice(pmSlice == 0) = NaN;
        
        % ------------------ Apply Requested Rotation ------------------
        % rotations(ori) times 90° CCW. Negative rotates clockwise.
        kRot = rotations(ori);
        if kRot ~= 0
            sourceSlice = rot90(sourceSlice, kRot);
            pmSlice     = rot90(pmSlice, kRot);
        end
        
        % ------------------ Plot the Slices ------------------
        subplot(opts.row, opts.col, jj);
        imagesc(sourceSlice);
        colormap(gray); freezeColors;
        axis off; axis image;

        subplot(opts.row, opts.col, jj+1);
        imagesc(pmSlice, opts.scale);
        colormap(map); freezeColors;
        axis off; axis image;
    end
    
    % Save figure for this orientation
    outBase = fullfile(foldername, [filename, '_ori', num2str(ori)]);
    f = gcf;
    
    % Save PNG
    exportgraphics(f, [outBase, '.png'], ...
        'Resolution', 600, 'BackgroundColor','black');
    % Save EPS
    exportgraphics(f, [outBase, '.eps'], ...
        'Resolution', 600, 'BackgroundColor','black');

    %close(f);
end

end % plotMap

% =============== Helper: getSlice ===============
function slice2D = getSlice(vol3D, mask3D, sliceIndex, orientation)
% Extract a 2D slice (plus any extra 4th-dim if present) according to:
%   ori=1 (transverse): vol3D(:,:,sliceIndex)
%   ori=2 (coronal)   : vol3D(:,sliceIndex,:)
%   ori=3 (sagittal)  : vol3D(sliceIndex,:,:)

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
