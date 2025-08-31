function [] = exportCurrentFig(path,filename)
%exports the current figure window in .png and .svg 
exportgraphics(gcf,[fullfile(path,filename),'.png'],'Resolution',600,'BackgroundColor','black');
print(gcf, [fullfile(path,filename),'.svg'], '-dsvg')
end