function [] = plotMap(sourceImg, mask, paramMap, map, filename, foldername, rotations, opts)
% plotMap: Plots parameter map next to source image in transverse, coronal, and sagittal views.
% Adds a dedicated colour‑bar column to the right (no change to existing subplot logic).
%
% Inputs and optional fields are identical to the original version.

% ------------------ Defaults ------------------
if ~isfield(opts,'scale'), opts.scale = [-5 5]; end
if ~isfield(opts,'row'),   opts.row   = 5;     end
if ~isfield(opts,'col'),   opts.col   = 6;     end   % TOTAL columns (image grid)

tight = 0;
if nargin<6 || isempty(foldername), foldername = pwd; end
if nargin<7 || isempty(rotations)  || numel(rotations)<3, rotations = [0 0 0]; end
filename   = char(filename);
foldername = char(foldername);

% ------------------ Prepare Inputs ------------------
sourceImg = double(sourceImg);
mask      = double(mask);  mask(mask > 0) = 1;
paramMap  = double(paramMap);

% ------------------ Loop through Orientations ------------------
for ori = 1:3
    % --- Compute non‑zero region in this orientation ---
    switch ori
        case 1, maskSum = squeeze(sum(sum(mask,1),2));   % Z
        case 2, maskSum = squeeze(sum(sum(mask,1),3));   % Y
        case 3, maskSum = squeeze(sum(sum(mask,2),3));   % X
    end
    nz = find(maskSum>0);
    if isempty(nz), warning('No data in orientation %d',ori); continue; end

    % focus on central 70 %% of mask slab
    rangeWidth = nz(end) - nz(1);            % scalar width of mask span
    slabStart  = round(nz(1)  + 0.15 * rangeWidth);
    slabEnd    = round(nz(end) - 0.15 * rangeWidth);

    % number of slice pairs (source + param)
    nSlices = max(1, floor((opts.row * opts.col) / 2));
    sliceIndices = round(linspace(slabStart, slabEnd, nSlices));

    % --- Create figure ---
    if tight
        figure('Color','k','Position',[10 100 100*opts.col 100*opts.row]);
    else
        figure('Color','k','Units','normalized','Position',[0.05 0.1 0.20 0.4]);   % widened by ~5 %
    end

    % ---- Plot each slice pair -----------------------------------------
    for s = 1:nSlices
        jj = (s-1)*2 + 1;          % left subplot index of pair
        idx = sliceIndices(s);

        srcSlice = getSlice(sourceImg, mask, idx, ori);
        prmSlice = getSlice(paramMap , mask, idx, ori);
        prmSlice = nanmean(prmSlice,4); prmSlice(prmSlice==0)=NaN;

        if rotations(ori)~=0
            srcSlice = rot90(srcSlice,rotations(ori));
            prmSlice = rot90(prmSlice,rotations(ori));
        end

        subplot(opts.row, opts.col, jj);
        imagesc(srcSlice); colormap(gray); freezeColors; axis off image;

        subplot(opts.row, opts.col, jj+1);
        imagesc(prmSlice, opts.scale); colormap(map); freezeColors; axis off image;
    end

    % ---- ADD EXTERNAL COLOUR‑BAR (does not affect subplot grid) -------
    cbAx = axes('Position',[0.88 0.1 0.02 0.8],'Visible','off');   % right side of figure
    colormap(cbAx, map); caxis(cbAx, opts.scale);
    cb = colorbar(cbAx, 'eastoutside');
    cb.Position = [0.91 0.1 0.02 0.8];   % fully outside image columns
    cb.Color = [1 1 1]; cb.Label.String = 'Parameter'; cb.Label.Color = [1 1 1];
    cb.Ticks = linspace(opts.scale(1), opts.scale(2), 5);

    % ---- Save output ---------------------------------------------------
    outBase = fullfile(foldername, sprintf('%s_ori%d',filename,ori));
    exportgraphics(gcf,[outBase '.png'],'Resolution',600,'BackgroundColor','black');
    exportgraphics(gcf,[outBase '.eps'],'Resolution',600,'BackgroundColor','black');
end
end

% ------------------ Helper: getSlice (unchanged) ------------------
function slice2D = getSlice(vol3D, mask3D, sliceIndex, orientation)
volSize = size(vol3D);
fallback = [1 1];

switch orientation
    case 1
        if sliceIndex>volSize(3), slice2D = nan(fallback); return; end
        slice = vol3D(:,:,sliceIndex,:);
        msk   = mask3D(:,:,sliceIndex);
    case 2
        if sliceIndex>volSize(2), slice2D = nan(fallback); return; end
        slice = vol3D(:,sliceIndex,:,:);
        msk   = mask3D(:,sliceIndex,:);
    case 3
        if sliceIndex>volSize(1), slice2D = nan(fallback); return; end
        slice = vol3D(sliceIndex,:,:,:);
        msk   = mask3D(sliceIndex,:,:);
    otherwise, slice2D = nan(fallback); return;
end

slice = squeeze(slice); msk = squeeze(msk);
if ndims(slice)==2
    slice2D = slice .* msk;
elseif ndims(slice)==3
    slice2D = bsxfun(@times, slice, msk);
else
    slice2D = nan(fallback);
end
end
