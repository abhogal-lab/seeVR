% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <ROIstd: calculates the standard deviation in a region of interest >
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
% ROImean calculates the mean value of the 'data' within an ROI specified
% by the mask
%
% data: 2D/3D data
%
% mask: corresponding 2D/3D mask
%
% val: average value within ROI
function [val] = ROIstd(data,mask)
mask = logical(mask);
data(isnan(data)) = 0; data(isinf(data)) = 0;
temp = data(:).*mask(:);
temp(temp == 0) = [];
val = std(temp);
end
