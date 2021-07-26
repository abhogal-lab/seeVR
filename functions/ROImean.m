%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function val = ROImean(data,mask)
%written by Alex Bhogal a.bhogal@umcutrecht.nl
%calculates the mean value in an ROI defined my the input mask

data(isnan(data)) = 0; data(isinf(data)) = 0;
temp = data(:).*mask(:);
temp(temp == 0) = [];
val = mean(temp);
end
