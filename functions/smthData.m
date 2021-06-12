function [sdata] = smthData(data, mask, opts)
%written by Alex Bhogal a.bhogal@umcutrecht.nl based on partial code from Jeroen
%Siero
global opts;
if isfield(opts,'spatialdim'); else; opts.spatialdim = 2; end 
if isfield(opts,'voxelsize'); else; opts.voxelsize = opts.headers.ts.dime.pixdim(2:4); end;
if isfield(opts,'filtWidth'); else; opts.filtWidth = 5; end %originally 7
if isfield(opts,'FWHM'); else; opts.FWHM = 4; end

sigma=opts.FWHM/2.355;
inplanevoxelsize=opts.voxelsize(1);
data_smooth=zeros(size(data));

filtSigma = sigma/inplanevoxelsize;
imageFilter=fspecial('gaussian',opts.filtWidth,filtSigma);

ndim = ndims(data);
%create a dilation mask
SE = strel('disk', 3);
dmask = imdilate(mask, SE);
delmask = double(dmask - mask); %cuts out center
ddata = zeros(size(data));

switch ndim
    case 3
        disp('Smoothing 3D image data' )
    data(isnan(data)) = 0;
    ddata = imdilate(data,SE);% dilate dataset to pad
    
    ddata = ddata.*delmask;
    ddata = ddata + data.*mask;
   
    switch opts.spatialdim
    case 2
        disp('Performing 2D spatial smoothing on image data' )
        for s=1:size(data, 3)
            if ~isempty(isnan(data))
                
                data_smooth(:,:,s)= nanconv(data(:,:,s),imageFilter, 'nanout');
            else
                data_smooth(:,:,s)=imgaussfilt(data(:,:,s),sigma/inplanevoxelsize, 'FilterDomain', 'spatial'); % here also fitler width of 7 is used internally
            end
        end
    case 3
        disp('Applying dilation to image to mitigate edge effects: check image edges' )
            disp('Performing 3D spatial smoothing on image data' )
            if isnan(data)
            
            ddata(isnan(ddata)) = 0;
            data_smooth=imgaussfilt3(ddata,sigma./opts.voxelsize, 'FilterDomain', 'spatial');
        else
            data_smooth=imgaussfilt3(ddata,sigma./opts.voxelsize, 'FilterDomain', 'spatial');
            
        end
    end
    case 4
        disp('Smoothing 4D timeseries data' )
    data(isnan(data)) = 0;
    
    for ii = 1:size(ddata,4)
        ddata(:,:,:,ii) = imdilate(data(:,:,:,ii),SE);% dilate dataset to pad
    end
    ddata = ddata.*delmask;
    ddata = ddata + data.*mask;
    
    switch opts.spatialdim
    case 2
        disp('Performing 2D spatial smoothing on timeseries data' )
        for t=1:size(data,4)
            for s=1:size(data, 3)
                data_smooth(:,:,s,t)= nanconv(data(:,:,s,t),imageFilter, 'nanout');
            end  
        end
        case 3
            disp('Applying dilation to image to mitigate edge effects: check image edges' )
            disp('Performing 3D spatial smoothing on timeseries data' )
        if isnan(data)
            ddata(isnan(ddata)) = 0;
            for t=1:size(data,4)
                data_smooth(:,:,:,t)=imgaussfilt3(ddata(:,:,:,t),sigma./opts.voxelsize, 'FilterDomain', 'spatial');
            end
        else
            %dilate image to mitigate edge effects
            for t=1:size(data,4)
                data_smooth(:,:,:,t)=imgaussfilt3(ddata(:,:,:,t),sigma./opts.voxelsize, 'FilterDomain', 'spatial');
            end
        end
    end
end
sdata = double(mask.*data_smooth);
opts.smoothmap = 1;
end
