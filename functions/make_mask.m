% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <grabTimeseres: extract voxel timeseries data and associated coordinates >
% 
% This program is free software: you can redistribute it and/or modify
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
% *************************************************************************
% function to generate a slice-by-slice brain mask based on user input
%
% input: input image to mask
%
% mask: resultant binary mask

function [mask] = make_mask(input)

slice = size(input,3);

for ii=1:slice
tmp = input(:,:,ii);
figure(1); imagesc(tmp); colormap gray 
h = impoly(); 
mask(:,:,slice) = createMask(h);

close(1)
end

disp('Finished creating mask');
end