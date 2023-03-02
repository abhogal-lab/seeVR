% Copyright (C) Alex A. Bhogal, 2021, University Medical Center Utrecht,
% a.bhogal@umcutrecht.nl
% <demeanData: subtracts the mean value of each voxel timeseries >
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
% This function calculates the mean value at each voxel and subtracts it
% from every element in the voxel timeseries. This can be used for 2-4D
% data.
%
% data: input timeseries data
%
% mask: binary mask defining region of interest
%
% dmData: a de-meaned version of the input data
function [dmData] = demeanData(data, mask)

warning('off')
global opts;
tf = class(data);
if isfield(opts,'verbose'); else; opts.verbose = 0; end %turn on/off select command output

switch ndims(data)
    
    case 2
        %if opts.verbose; disp('de-meaning timeseries vector'); end
        m_ts = mean(data);
        dmData = data-m_ts;
    case 3
        if opts.verbose; disp('de-meaning 2D timeseries data'); end
        [x,y,dyn] = size(data);
        [voxel_ts,coordinates] = grabTimeseries(data,mask);
        m_ts = nanmean(voxel_ts,2);
        voxel_ts = voxel_ts-m_ts;
        dmData = zeros(size(data));
        dmData = reshape(dmData,[x*y dyn]);
        dmData(coordinates,:) = voxel_ts;
        dmData = reshape(dmData,size(data));
    case 4
        %if opts.verbose; disp('de-meaning 3D timeseries data'); end
        [x,y,z,dyn] = size(data);
        [voxel_ts,coordinates] = grabTimeseries(data,mask);
        m_ts = nanmean(voxel_ts,2);
        voxel_ts = voxel_ts-m_ts;
        dmData = zeros(size(data));
        dmData = reshape(dmData,[x*y*z dyn]);
        dmData(coordinates,:) = voxel_ts;
        dmData = reshape(dmData,size(data));
end
dmData = cast(dmData,tf);
end