% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <view_ts: a function to view timeseries data >
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
% This function that displays timeseries data
%
% data: input 4D timeseries
% slice: which slice of the 3D volume to show
% scale: range of the data to be shown
function [] = view_ts(data,slice,scale)


figure;
for ii=1:size(data,4)
    imagesc(imrotate(data(:,:,slice,ii),90), scale)
    colormap(jet)
    pause(0.01)
end

end