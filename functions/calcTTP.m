function [ttp_map] = calculateTTPfromtau(data, probe, mask, opts)

global opts
input_probe = rescale(probe);

coordinates = find(logical(mask(:)));

% Smooth the CO2 probe data
window_size = 5;
smoothed_probe = input_probe;

%calculate highest positive signal change (corresponding to first
%onset)
dt_smoothed_probe = rescale(gradient(smoothed_probe), -1, 1);
figure;
subplot(2,1,1); plot(rescale(smoothed_probe)); hold on; plot(dt_smoothed_probe)
[PKS,LOCS] = findpeaks(dt_smoothed_probe,'MinPeakHeight',0.4,'NPeaks',3)
subplot(2,1,2); plot(dt_smoothed_probe,'b'); hold on; scatter(LOCS(1),PKS(1));
%first higher peak will correspond to rising edge
on_index = LOCS(1) %start search around here
%do the same for the negative signal
[PKS,LOCS] = findpeaks(-dt_smoothed_probe,'MinPeakHeight',0.4,'NPeaks',3)
hold on;  plot(-dt_smoothed_probe, 'g'); hold on; scatter(LOCS(1),PKS(1));
off_index = LOCS(1)

threshold = max(dt_smoothed_probe) * 0.1; % 10% of maximum gradient

% find the value at which the signal crosses the threshold

for ii=(on_index-15):on_index
    if dt_smoothed_probe(1,ii) > threshold
        onset = ii-1
        break
    else
        continue
    end
end

%peak search window
window = [onset:1:off_index+10];

%crop data
cr_data = data(:,window);
%interpolate data by factor 10
parfor ii = 1:length(data)
    iCr_data(ii,:) = interp(cr_data(ii,:),10);
end

ttp_values = zeros(length(data),1);
for ii=1:length(cr_data)
    [~,I] = max(iCr_data(ii,:));
    ttp_values(ii) = (I/10)*opts.TR;
end

ttp_map = zeros(1, numel(mask));
ttp_map(coordinates) = ttp_values;
ttp_map = reshape(ttp_map, size(mask));

end