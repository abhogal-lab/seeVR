%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Carpet Plot Edge Detection
% Created By: Brad Fitzgerald - Purdue University, fitzge45@purdue.edu
% Last Edited On: 12 Nov 2020

% Created for publication with corresponding manuscript, "Using carpet 
% plots to analyze transit times of low frequency oscillations in resting
% state fMRI" - submitted to Scientific Reports


% This script is designed to analyze some input matrix which represents a
% carpet plot of fMRI data (where voxels are ordered along the y-axis and
% time-series volumes are ordered along the x-axis) and detect edges formed
% by low frequency oscillations in the fMRI BOLD signal. The script was
% designed to analyze carpet plots which have previously had the voxels
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
clear
close all;

%Input carpet plot matrix
im1 = input_img;
TR = 0.72;

%Choose number of lines to look for
num_lines = 25;

%Choose to calc on rising (black to white) or falling (white to black)
%edges
%-1 = rising edge, 1 = falling edge
edge = -1;

%Set contrast thereshold for limiting accepted edges
contrast_threshold = 0.2;

%Set neuronal signal threshold
neur_thresh = 0.15;

%% Smooth image

scale_factor = 1; %Used to scale computation of edge angles if desired
height = size(im1,1);

orig_avg = mean(im1);

sort_avg = sort(orig_avg, 'descend');
neur_thresh_val = sort_avg(floor(size(sort_avg,2)*neur_thresh));

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
old_im1 = im1;
im1 = smoothed_data;

%% Average Time Series

%Find and plot average time series
avg = mean(im1(:,:));
% figure
% plot(avg);
% xlabel('fMRI volume');
% ylabel('Average across voxels');
% title('Average time series');

smoothed_avg = avg;

%% Derivative
%Apply derivative filter
h = 1/2*[1 0 -1];
filtsize = floor(size(h,2) / 2);
derivatived_avg = edge * filter2(h, smoothed_avg);
%Adjust the ends of the derivatived data since the front and back will be
%skewed
for i=1:filtsize
    derivatived_avg(filtsize-(i-1)) = derivatived_avg(filtsize-(i-2));
    derivatived_avg(size(derivatived_avg,2)-filtsize+(i)) = derivatived_avg(size(derivatived_avg,2)-filtsize+(i-1));
end

% figure
% plot(derivatived_avg)
% title('Derivatived Data')

%% Find Peaks
max_ind = zeros(1,num_lines);
temp = derivatived_avg;

%Here we find the n=num_lines highest peaks of the derivative data,
%representing locations where we want to draw a line
for i=1:num_lines
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
        num_lines = i
        break;   
    end
end

neur_assignment = zeros(1, num_lines);

max_ind = sort(max_ind(max_ind>0), 'ascend');
%Plot derivative over data to illustrate where the algorithm is deciding to
%draw lines
% figure
% imagesc(im1);
% hold on;
% colormap(gray);
% caxis([-1 1]);
% plot(derivatived_avg*20000, 'r', 'LineWidth', 2);
% hold off;

%% Draw Lines
%Initialize space to store line data
linfit = zeros(num_lines, 2);
linfitline = zeros(num_lines, height-1+1);

%Apply derivative filter to all data
h = 1/8*[1 0 -1; 2 0 -2; 1 0 -1];
filtsize = floor(size(h,2) / 2);
derivatived_data = edge * filter2(h, im1);
for b=1:filtsize
    derivatived_data(:,filtsize-(b-1)) = derivatived_data(:,filtsize-(b-2));
    derivatived_data(:,size(derivatived_data,2)-filtsize+(b)) = derivatived_data(:,size(derivatived_data,2)-filtsize+(b-1));
end

temp = derivatived_avg;
avg_contrast = zeros(1, num_lines);
for i=1:num_lines
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
    %We can display the considered data if we want, but it's commented out
    %as seen below
%     figure
%     imshow(data);
%     colormap(gray);    
    
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
%     figure
    %plot(flip(max_ind2), datayax);
%     scatter(flip(max_ind2), datayax);
%     hold on;
    linfit(i,:) = polyfit(datayax, max_ind2, 1);
    linfitline(i,:) = linfit(i,1)*datayax + linfit(i,2);
%     plot(flip(linfitline), datayax);
%     hold off;
%     title('Finding Best Fit Line')
    
end

%Calculate slopes and angles
slopes = -1./linfit(:,1);
angles = atan(slopes * scale_factor)*180/pi;
times = ones(size(slopes)) * size(im1,1) ./ slopes;

%% Main Fig 1: Create figure with all edges and transit times for this subject

figure('Position', [50 50 600 400])
h1 = subplot(2,1,1);
imagesc(old_im1);
caxis([-1 1]);
colormap(gray);
hold on;
plot_mat = [];
for i=1:num_lines
    if avg_contrast(i) > contrast_threshold
        if slopes(i) < 0
            plot((linfitline(i,:)-mean(linfitline(i,:))+max_ind(i)), datayax+1, 'g', 'LineWidth', 2);
        else
            plot((linfitline(i,:)-mean(linfitline(i,:))+max_ind(i)), datayax+1, 'r', 'LineWidth', 2);
        end
        plot_mat = [plot_mat; [times(i)*TR max_ind(i)*TR]];
    end
end
yticklabels([]);
xticklabels([]);
ylabel('Voxels', 'fontweight', 'bold', 'FontSize', 12);
title(strcat('Subject', {' '}, string(sub)), 'FontSize', 14);

hold off;



h1 = subplot(2,1,2);
h1.Position = h1.Position + [0 0.185 0 -0.075];
plot(plot_mat(:,2), plot_mat(:,1), '-ok', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
xlim([1 360]);
xticks([50 100 150 200 250 300 350]);
ylim([-4 10]);
yticks([-4 0 4 8]);
ax = gca;
ax.FontSize = 12;
grid on;
yline(0, '--');
xlabel('Time (s)', 'fontweight', 'bold', 'FontSize', 12);
ylabel('Transit time (s)', 'fontweight', 'bold', 'FontSize', 12);

%% Main Fig 2: Create figure with neuro assigned edges and average BOLD signal

figure('Position', [50 50 800 500])

h1 = subplot(2,1,1);
h1.Position = h1.Position + [0 0.15 0 -0.15];
plot(orig_avg, 'k');
hold on;
yline(neur_thresh_val, 'r--');
ylabel('Avg. BOLD signal', 'Fontweight', 'bold');
xticklabels([]);
title(strcat('Carpetplot edges at top', {' '}, string(neur_thresh*100), '% of average BOLD time series'));

L = plot(nan, nan, 'r--');
hLegend = legend(L, {strcat(string((1-neur_thresh)*100), '% Threshold')});
hLegend.Location = 'southeast';

hold off;

h2 = subplot(2,1,2)
h2.Position = h2.Position + [0 0 0 0.27];

imagesc(old_im1);
colormap(gray);
caxis([-1 1]);
hold on;

for i=1:num_lines
    if avg_contrast(i) > 0.2 && neur_assignment(i) == 1
        plot((linfitline(i,:)-mean(linfitline(i,:))+max_ind(i)), datayax+1, 'b', 'LineWidth', 2);
    elseif avg_contrast(i) > 0.2 && neur_assignment(i) == 2
        plot((linfitline(i,:)-mean(linfitline(i,:))+max_ind(i)), datayax+1, 'r', 'LineWidth', 2);
    else
        neur_assignment(i) = 0;
    end

end
hold off;

yticklabels([]);
xlabel('Time (s)', 'fontweight', 'bold');
xticks([50:50:500]);
xticklabels({'36', '72', '108', '144', '180', '216', '252', '288', '324', '360'});
ylabel('Voxels', 'fontweight', 'bold');
hold on;
L2(1) = plot(nan, nan, 'b-');
L2(2) = plot(nan, nan, 'r-');
hLegend2 = legend(L2, {strcat('Lower', {' '}, string(100 - neur_thresh*100), '%'), strcat('Top', {' '}, string(neur_thresh*100), '%')});
hLegend2.Location = 'southeast';
hold off;
