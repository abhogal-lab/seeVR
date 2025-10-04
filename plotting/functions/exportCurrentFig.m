function exportCurrentFig(path, filename, fig)
% exports a specific figure window to .png and .svg
if nargin < 3 || ~ishandle(fig) || ~isgraphics(fig)
    fig = gcf;  % fallback
end
set(fig,'InvertHardCopy','off','WindowState','normal');  % reliable sizing
drawnow;  % ensure the figure is fully rendered
%exportgraphics(fig, [fullfile(path,filename) '.png'], 'Resolution',600, 'BackgroundColor','black');
print(fig, [fullfile(path,filename) '.png'], '-dpng', '-r600');   % PNG
print(fig, [fullfile(path,filename) '.svg'], '-dsvg');
end
