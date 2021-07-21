%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [mask] = make_mask(input)
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