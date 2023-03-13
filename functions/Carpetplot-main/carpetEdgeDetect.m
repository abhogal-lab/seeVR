% Modifications to original code are Copyright (C) Alex A. Bhogal, 2022, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl. Please see the original license file provided by
% Tong et. al that is included with this code.
%
% <carpetEdgeDetect: performs shag carpet plot analysis and outputs transit maps >
%
% This modified program is free software: you can redistribute it and/or modify
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
function [maps] = carpetEdgeDetect(data,sort_map, mask, opts)

global opts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Carpet Plot Edge Detection
% Created By: Brad Fitzgerald - Purdue University, fitzge45@purdue.edu
% Last Edited On: 12 Nov 2020

% Created for publication with corresponding manuscript, "Using carpet
% plots to analyze transit times of low frequency oscillations in resting
% state fMRI" - submitted to Scientific Reports

% Updated for integration with seeVR toolbox (www.seeVR.nl)
% Modified by: Alex Bhogal - UMC Utrecht, a.bhogal@umcutrecht.nl

% This script is designed to analyze some input matrix which represents a
% carpet plot of fMRI data (where voxels are ordered along the y-axis and
% time-series volumes are ordered along the x-axis) and detect edges formed
% by low frequency oscillations in the fMRI BOLD signal. The script was
% designed to anFalyze carpet plots which have previously had the voxels
% (rows) specifically ordered based on arrival time of the low frequency
% oscillation pattern in each voxel - this helps to form more clearly
% defined "edges" in the carpet plot. This was the intended use of the
%script, though edge detection could technically be performed without this
%ordering (but results may be poor). The script also was built for carpet
%plots where the data has been demeaned and scaled down according to
%standard deviation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Figures:

% Main Figure 1: a plot which shows the detected edges overlayed on the
% carpet plot with computed edge widths (transit time) plotted below

% Main Figure 2: a plotwhich shows the detected edges separated based on
% whether the edges correspond to the top 15% (or some redefined threshold)
% of the voxel-averaged time series.

% Additional figures are written to help analyze the computation process,
% but are currently commented out - they can be uncommented if you wish to
% use them.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Inputs:

% im1 = input carpet plot matrix

% TR = TR from fMRI imaging sequence

% num_lines = number of edges to try to find (will stop at highest number
%                       detectable if unable to reach this number)

% edge = decides whether you want to detect edges of increasing signal
%                       intensity or decreasing signal intesity

% contrast_threshold = decides how strong the change in signal intensity
%must be for an edge to be accepted

% neur_thresh = threshold for dividing edges for Main Fig 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs
%clear
%close all;

%Input carpet plot matrix
%im1 = input_img;
%TR = 0.72;

%Choose number of lines to look for
%num_lines = 25;

%Choose to calc on rising (black to white) or falling (white to black)
%edges
%-1 = rising edge, 1 = falling edge
%edge = -1;

%Set contrast thereshold for limiting accepted edges
%contrast_threshold = 0.2;

%Set neuronal signal threshold
%neur_thresh = 0.15;
%set defaults

if isfield(opts.carpet,'num_lines'); else; opts.carpet.num_lines = 25; end
if isfield(opts.carpet,'edge'); else; opts.carpet.edge = -1; end
if isfield(opts.carpet,'contrast_threshold'); else; opts.carpet.contrast_threshold = .2; end
if isfield(opts.carpet,'neur_thresh'); else; opts.carpet.neur_thresh = 0.15; end
if isfield(opts.carpet,'rescale'); else; opts.carpet.rescale = 1; end
if isfield(opts.carpet,'knots'); else; opts.carpet.knots = 3; end
if isfield(opts.carpet,'fittype'); else; opts.carpet.fittype = 'spline'; end                    %can be 'spline' or 'linear'
if isfield(opts.carpet,'old_smooth'); else; opts.carpet.old_smooth = 0'; end                    %added smoothing option
if isfield(opts.carpet,'invert'); else; opts.carpet.invert = 0'; end                            %inversion can help with certain datasets (i.e. Gd)
if isfield(opts.carpet,'peak'); else; opts.carpet.peak = 0'; end                                %calculate based on inflections (i.e. Gd)
if isfield(opts.carpet,'interp_factor'); else; opts.carpet.interp_factor = 1; end                      %factor by which to temporally interpolate data.


%setup savedir
opts.carpetdir = fullfile(opts.resultsdir,'carpet'); mkdir(opts.carpetdir);

% delete any existing transit maps
cd(opts.carpetdir);
delete *transit*
%% sort input data

sort_map = sort_map(:);
data = demeanData(data,mask);
[voxels,coords] = grabTimeseries(data , mask);
ivoxels = zeros([length(coords), size(voxels,2)*opts.carpet.interp_factor]);

%interpolate data
parfor ii = 1:length(coords)
    ivoxels(ii,:) = interp(voxels(ii,:),opts.carpet.interp_factor);
end
voxels = ivoxels; clear ivoxels;

dvoxels = zeros([size(voxels,1) size(voxels,2)-1]);

if opts.carpet.peak
    parfor ii=1:size(voxels,1)
        dvoxels(ii,:) = diff(voxels(ii,:));
    end
    dvoxels(:,end+1) = dvoxels(:,end);
    voxels = dvoxels; clear dvoxels;
end

if opts.carpet.rescale
    if opts.carpet.invert
        parfor ii=1:size(voxels,1)
            voxels(ii,:) = abs(rescale(voxels(ii,:))-1);
        end
    else
        parfor ii=1:size(voxels,1)
            voxels(ii,:) = rescale(voxels(ii,:));
        end
    end
else
    if opts.carpet.invert
        parfor ii=1:size(voxels,1)
            voxels(ii,:) = abs(voxels(ii,:)-max(voxels(ii,:)));
        end
    end
end


%sorted index at all non-zero coordinates
[~,I] = sort(sort_map(coords), 'descend');
im1 = voxels(I,:);

figure;
subplot(1,2,1);
imagesc(imgaussfilt(voxels)); colormap(gray); title('unsorted')
subplot(1,2,2);
imagesc(imgaussfilt(im1)); colormap(gray); title('sorted')
saveas(gcf,[opts.carpetdir,'sorting.fig']);

%generate a coordinate mask
carpet_mask = zeros([1 numel(mask)]);
carpet_mask(coords) = 1;
carpet_mask = reshape(carpet_mask,size(mask));
if opts.niiwrite
    cd(opts.carpetdir)
    niftiwrite(cast(carpet_mask, opts.maskDatatype),'carpetMask',opts.info.mask);
else
    saveImageData(carpet_mask, opts.headers.mask, opts.carpetdir, 'carpetMask.nii.gz', 64);
end
maps.carpetMask = carpet_mask;

%% Smooth image
% if opts.GD
% im1(im1>0) = 1;
% end

xdata = opts.TR:opts.TR:opts.TR*size(im1,2);
scale_factor = 1; %Used to scale computation of edge angles if desired
height = size(im1,1);

orig_avg = mean(im1);

sort_avg = sort(orig_avg, 'descend');
neur_thresh_val = sort_avg(floor(size(sort_avg,2)*opts.carpet.neur_thresh));

if opts.carpet.old_smooth
    %Create filter h for image blurring - meant to smooth out image for
    %clearer edges
    h = 1/12*[1 1.5 1;
        1.5 2 1.5;
        1 1.5 1];
    h = conv2(h,h);
    
    %Filter data for several repititions
    smoothed_data = filter2(h, im1);
    for i=1:5
        smoothed_data = filter2(h, smoothed_data);
    end
else
    opts.filter = 'gaussian';
    opts.FMWH = [12 12];
    opts.spatialdim = 2;
    im_mask = ones(size(im1));
    [smoothed_data] = filterData(im1,im1,im_mask,opts);
end

old_im1 = im1;
im1 = smoothed_data;

%figure; colormap(gray); title('smoothing'); subplot(2,1,1); imagesc(old_im1); subplot(2,1,2); imagesc(smoothed_data)
%% Average Time Series

%Find and plot average time series
avg = mean(im1(:,:));
smoothed_avg = avg;

%% Derivative
%Apply derivative filter
h = 1/2*[1 0 -1];
filtsize = floor(size(h,2) / 2);
derivatived_avg = opts.carpet.edge * filter2(h, smoothed_avg);
%Adjust the ends of the derivatived data since the front and back will be
%skewed
for i=1:filtsize
    derivatived_avg(filtsize-(i-1)) = derivatived_avg(filtsize-(i-2));
    derivatived_avg(size(derivatived_avg,2)-filtsize+(i)) = derivatived_avg(size(derivatived_avg,2)-filtsize+(i-1));
end

%% Find Peaks
max_ind = zeros(1,opts.carpet.num_lines);
temp = derivatived_avg;

%Here we find the n=num_lines highest peaks of the derivative data,
%representing locations where we want to draw a line
for i=1:opts.carpet.num_lines
    is_line_on_back_edge = 1;
    is_line_on_front_edge = 1;
    while is_line_on_back_edge == 1 || is_line_on_front_edge == 1
        [~, max_ind(i)] = max(temp);
        a=0;
        while (max_ind(i) + a <= size(im1,2)) && (temp(max_ind(i) + a) > 0)
            temp(max_ind(i) + a) = 0;
            a= a+1;
        end
        if max_ind(i) + a <= size(im1,2)
            is_line_on_back_edge = 0;
        else
            is_line_on_back_edge = 1;
        end
        a=-1;
        while (max_ind(i) + a >= 1) && (temp(max_ind(i) + a) > 0)
            temp(max_ind(i) + a) = 0;
            a= a-1;
        end
        if max_ind(i) + a >= 1
            is_line_on_front_edge = 0;
        else
            is_line_on_front_edge = 1;
        end
    end
    if nnz(temp<=0) == length(temp)
        disp('Not enough locations for a line!');
        opts.carpet.num_lines = i;
        break;
    end
end

neur_assignment = zeros(1, opts.carpet.num_lines);

max_ind = sort(max_ind(max_ind>0), 'ascend');

%% Draw Lines
%Initialize space to store line data
linfit = zeros(opts.carpet.num_lines, 2);
linfitline = zeros(opts.carpet.num_lines, height-1+1);

%Apply derivative filter to all data
h = 1/8*[1 0 -1; 2 0 -2; 1 0 -1];
filtsize = floor(size(h,2) / 2);
derivatived_data = opts.carpet.edge * filter2(h, im1);
for b=1:filtsize
    derivatived_data(:,filtsize-(b-1)) = derivatived_data(:,filtsize-(b-2));
    derivatived_data(:,size(derivatived_data,2)-filtsize+(b)) = derivatived_data(:,size(derivatived_data,2)-filtsize+(b-1));
end

temp = derivatived_avg;
avg_contrast = zeros(1, opts.carpet.num_lines);
for i=1:opts.carpet.num_lines
    location = max_ind(i);
    fin(i)=0;
    %For now we are considering the range around the peak where the
    %derivative is greater than 0, plus adding an extra 2 points on either
    %side
    while (max_ind(i) + fin(i) <= size(im1,2)) && (temp(max_ind(i) + fin(i)) > 0)
        temp(max_ind(i) + fin(i)) = 0;
        fin(i)= fin(i)+1;
    end
    st(i)=-1;
    while (max_ind(i) + st(i) >= 1) && (temp(max_ind(i) + st(i)) > 0)
        temp(max_ind(i) + st(i)) = 0;
        st(i)= st(i)-1;
    end
    
    %Compute averaged contrast
    cont_st = smoothed_avg(st(i)+max_ind(i));
    cont_end = smoothed_avg(fin(i)+max_ind(i));
    avg_contrast(i) = abs(cont_end - cont_st);
    
    %Assign edges for Main Fig 2
    if orig_avg((max_ind(i)+fin(i))) > neur_thresh_val
        neur_assignment(i) = 2;
    else
        neur_assignment(i) = 1;
    end
    
    fin(i) = fin(i)+ 2;
    st(i) = st(i) - 2;
    
    st(i)=st(i)+max_ind(i);
    if st(i) < 1
        st(i) = 1;
    end
    fin(i)=fin(i)+max_ind(i);
    if fin(i) > size(im1,2)
        fin(i) = size(im1,2);
    end
    data = derivatived_data(1:height, st(i):fin(i));
    
    %Find horizontal point where derivative is maximized for each row
    max_ind2 = zeros(1, size(data, 1));
    for b=1:size(data,1)
        [~, max_ind2(b)] = max(data(b, :));
    end
    %Account for extreme values that will probably throw off algorithm
    max_ind2(length(max_ind2)) = max_ind2(length(max_ind2)-1);
    max_ind2(1) = max_ind2(2);
    
    %Calculate the best fit line for our derivative peak locations
    datayax = 1:size(max_ind2,2);
    
    switch opts.carpet.fittype
        case 'linear'
            linfit(i,:) = polyfit(datayax, max_ind2, 1);
            linfitline(i,:) = linfit(i,1)*datayax + linfit(i,2);
        case 'spline' %default
            slm = slmengine(datayax, max_ind2, 'knots', opts.carpet.knots);
            linfitline(i,:) = slmeval(datayax, slm);
    end
    
end

%Calculate slopes and angles
switch opts.carpet.fittype
    case 'linear'
        slopes = -1./linfit(:,1);
        angles = atan(slopes * scale_factor)*180/pi;
        times = ones(size(slopes)) * size(im1,1) ./ slopes;
    case 'spline' %default
        slopes = linfitline(:,1) - linfitline(:,end)
        times = slopes;
end

%% Main Fig 1: Create figure with all edges and transit times for this subject

figure('Position', [50 50 600 400])
h1 = subplot(2,1,1);
imagesc(old_im1);
caxis([-1 1]);
colormap(gray);
hold on;
plot_mat = [];
mm = 0;
for i=1:opts.carpet.num_lines
    if avg_contrast(i) > opts.carpet.contrast_threshold
        mm=mm+1;
        if slopes(i) < 0
            data_line(i,:) = (linfitline(i,:)-mean(linfitline(i,:))+max_ind(i));
            transit_line(mm,:) = data_line(i,:) - data_line(i,1);
            transit_line(mm,1) = transit_line(mm,2)/2;
            plot(data_line(i,:), datayax+1, 'g', 'LineWidth', 2);
        else
            data_line(i,:) = (linfitline(i,:)-mean(linfitline(i,:))+max_ind(i));
            transit_line(mm,:) = data_line(i,:) - data_line(i,end);
            transit_line(mm,1) = transit_line(mm,2)/2;
            plot(data_line(i,:), datayax+1, 'r', 'LineWidth', 2);
        end
        plot_mat = [plot_mat; [times(i)*opts.TR max_ind(i)*opts.TR]];
    end
end
yticklabels([]);
xticklabels([]);
ylabel('Voxels', 'fontweight', 'bold', 'FontSize', 12);

hold off;

h1 = subplot(2,1,2);
h1.Position = h1.Position + [0 0.185 0 -0.075];
plot(plot_mat(:,2), plot_mat(:,1)/opts.carpet.interp_factor, '-ok', 'MarkerFaceColor', 'k', 'MarkerSize', 3);

ax = gca;
ax.FontSize = 12;
grid on;
yline(0, '--');
xlabel('Time (s)', 'fontweight', 'bold', 'FontSize', 12);
ylabel('Transit time (s)', 'fontweight', 'bold', 'FontSize', 12);

if opts.carpet.edge > 0
    saveas(gcf,[opts.carpetdir,'times_falling.fig']);
else
    saveas(gcf,[opts.carpetdir,'times_rising.fig']);
end

%generate transit maps

if opts.niiwrite
    cd(opts.carpetdir)
    niftiwrite(cast(carpet_mask, opts.maskDatatype),'carpetMask',opts.info.mask);
else
    saveImageData(carpet_mask, opts.headers.mask, opts.carpetdir, 'carpetMask.nii.gz', 64);
end

for ii=1:mm
    transit_map = zeros([1 numel(mask)]);
    transit_map(coords(flip(I))) = opts.TR*flip(transit_line(ii,:));
    transit_map = reshape(transit_map,size(mask));
    if any(transit_line(ii,:)<0)
        transit_map(transit_map ~= 0) =  transit_map(transit_map ~= 0) + abs(min(transit_line(ii,:)));
    end
    if opts.niiwrite
        cd(opts.carpetdir)
        niftiwrite(cast(transit_map/opts.carpet.interp_factor,opts.mapDatatype),['transitMap_',int2str(ii)],opts.info.map);
    else
        saveImageData(transit_map/opts.carpet.interp_factor, opts.headers.map, opts.carpetdir, ['transitMap_',int2str(ii),'.nii.gz'], 64);
    end
    mapNr = ['transitMap_',int2str(ii)];
    eval(['maps.',mapNr,' = transit_map']);
end

%generate average transit map
mean_transit = mean(transit_line,1);
transit_map = zeros([1 numel(mask)]);
transit_map(coords(flip(I))) = opts.TR*flip(mean_transit);
transit_map = reshape(transit_map,size(mask));
if any(transit_line(ii,:)<0)
    transit_map(transit_map ~= 0) =  transit_map(transit_map ~= 0) + abs(min(transit_line(ii,:)));
end
if opts.niiwrite
    cd(opts.carpetdir)
    niftiwrite(cast(transit_map/opts.carpet.interp_factor, opts.mapDatatype),'meanTransitMap',opts.info.map);
else
    saveImageData(transit_map/opts.carpet.interp_factor, opts.headers.map, opts.carpetdir, 'meanTransitMap.nii.gz', 64);
end
maps.meanTransit = transit_map/opts.carpet.interp_factor;
%% Main Fig 2: Create figure with neuro assigned edges and average BOLD signal

figure('Position', [50 50 800 500])

h1 = subplot(2,1,1);
h1.Position = h1.Position + [0 0.15 0 -0.15];
plot(xdata, orig_avg, 'k');
hold on;
%plot(plot_mat(:,2), plot_mat(:,1), '*', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
xlim([0 xdata(end)]);

hold on;
yline(neur_thresh_val, 'r--');
ylabel('Avg. BOLD signal', 'Fontweight', 'bold');
xticklabels([]);
title(strcat('Carpetplot edges at top', {' '}, string(opts.carpet.neur_thresh*100), '% of average BOLD time series'));

L = plot(nan, nan, 'r--');
hLegend = legend(L, {strcat(string((1-opts.carpet.neur_thresh)*100), '% Threshold')});
hLegend.Location = 'southeast';

hold off;

h2 = subplot(2,1,2)
h2.Position = h2.Position + [0 0 0 0.27];

imagesc(old_im1);
colormap(gray);
caxis([-1 1]);
hold on;

for i=1:opts.carpet.num_lines
    if avg_contrast(i) > 0.2 && neur_assignment(i) == 1
        plot(data_line(i,:), datayax+1, 'b', 'LineWidth', 2);
    elseif avg_contrast(i) > 0.2 && neur_assignment(i) == 2
        plot(data_line(i,:), datayax+1, 'r', 'LineWidth', 2);
    else
        neur_assignment(i) = 0;
    end
    
end
hold off;

yticklabels([]);
xlabel('Time (s)', 'fontweight', 'bold');
ylabel('Voxels', 'fontweight', 'bold');
hold on;
L2(1) = plot(nan, nan, 'b-');
L2(2) = plot(nan, nan, 'r-');
hLegend2 = legend(L2, {strcat('Lower', {' '}, string(100 - opts.carpet.neur_thresh*100), '%'), strcat('Top', {' '}, string(opts.carpet.neur_thresh*100), '%')});
hLegend2.Location = 'southeast';
hold off;

if opts.carpet.edge > 0
    saveas(gcf,[opts.carpetdir,'lines_falling.fig']);
else
    saveas(gcf,[opts.carpetdir,'lines_rising.fig']);
end

%save maps
disp('saving maps in .mat file' )
carpet_maps = maps;
save([opts.carpetdir,'carpet_Maps.mat'], 'carpet_maps');

end
