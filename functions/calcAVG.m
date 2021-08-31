% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <calcAVG:  calculates the mean value of a vector between specified points>
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
% function to calculate the mean value of a vector between specified points
%
% data: 1D vector
%
% idx: user supplied index (optional)
% if no index is specified use can select start and end point themselves
function [idx, avgVal] = calcAVG(vector, idx)

if iscolumn(vector); else; vector = vector'; end
switch nargin
    case 2 %user supplied index
        avgVal = mean(vector(idx(1):dx(2),1))
    case 1 % user defined points
        figure; plot(vector)
        title('Select start and end point of epoch');
        [idx,~] = ginput(2); idx = round(idx);
        close;
        avgVal = mean(vector(idx(1):idx(2),1));
end
end