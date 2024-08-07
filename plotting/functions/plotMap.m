function [] = plotMap(sourceImg,mask,paramMap,map,opts)
% a.bhogal@umcutrecht.nl
% <plotMap: plots parameter map next to source image >
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
%
% sourceImg: source image (e.g. unprocessed mean BOLD image or anatomical image in the case that it has the same matrix size as the parameter map)
%
% mask: binary mask defining voxels of interest
%
% paramMap: parameter map with same marix size as sourceImg
%
% map: colormap to be used
%
% opts: options structure containing required variables for this specific
% function; i.e. opts.scale opts.row, opts.col, opts,niiwrite, opts.resultsdir
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(opts,'scale'); else; opts.scale = [-5 5]; end                   % data range
if isfield(opts,'row'); else; opts.row = 5; end                            % nr rows in plot
if isfield(opts,'col'); else; opts.col = 6; end                            % nr cols in plot
if isfield(opts,'step'); else; opts.step = 1; end                          % nr of steps through slices
if isfield(opts,'start'); else; opts.start = floor(size(mask,3))/3; end    % first slice number


sourceImg = double(sourceImg);
mask = double(mask);
paramMap = double(paramMap);

    figure;
    set(gcf, 'Units', 'pixels', 'Position', [10, 100, 1500, 900]);

    for jj=1:2:opts.row*opts.col-1  
    subplot(opts.row,opts.col,jj); imagesc(mask(:,:,opts.start + jj*opts.step).*sourceImg(:,:,opts.start + jj*opts.step)); 
    colormap(gray); freezeColors; set(gca,'visible','off'); set(gca,'xtick',[]);set(gcf,'color','k') 
        
    tmp = paramMap(:,:,opts.start + jj*opts.step);
    tmp = nanmean(tmp,4).*mask(:,:,opts.start + jj*opts.step); tmp(tmp == 0)= NaN;
    subplot(opts.row,opts.col,jj+1); imagesc(tmp,opts.scale); colormap(map); 
    freezeColors; set(gca,'visible','off'); set(gca,'xtick',[]);set(gcf,'color','k') 

    end
      
end


