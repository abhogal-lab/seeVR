function [mask] = make_mask(input)
%%written by Alex Bhogal a.bhogal@umcutrecht.nl
% function to generate a slice-by-slice brain mask

slice = size(input,3)

for ii=1:slice
tmp = input(:,:,ii);
figure(1); imagesc(tmp); colormap gray 
h = impoly(); 
mask(:,:,slice) = createMask(h);

close(1)
end

disp('Finished creating mask');
end