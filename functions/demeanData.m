%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht,
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes.

function [dmData] = demeanData(data, mask)
switch ndims(data)
    
    case 2
        disp('de-meaning timeseries vector')
        m_ts = mean(data);
        dmData = data-m_ts;
    case 3
        disp('de-meaning 2D timeseries data')
        [x,y,dyn] = size(data);
        [voxel_ts,coordinates] = grabTimeseries(data,mask);
        m_ts = nanmean(voxel_ts,2);
        voxel_ts = voxel_ts-m_ts;
        dmData = zeros(size(data));
        dmData = reshape(dmData,[x*y dyn]);
        dmData(coordinates,:) = voxel_ts;
        dmData = reshape(dmData,size(data));
    case 4
        disp('de-meaning 3D timeseries data')
        [x,y,z,dyn] = size(data);
        [voxel_ts,coordinates] = grabTimeseries(data,mask);
        m_ts = nanmean(voxel_ts,2);
        voxel_ts = voxel_ts-m_ts;
        dmData = zeros(size(data));
        dmData = reshape(dmData,[x*y*z dyn]);
        dmData(coordinates,:) = voxel_ts;
        dmData = reshape(dmData,size(data));

end

end