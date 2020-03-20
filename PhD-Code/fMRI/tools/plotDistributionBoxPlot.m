function plotDistributionBoxPlot(data,origin,label)

color = [255 0 0]./255;

f = figure;

boxplot(data);
h = findobj(gca,'tag','Median');
set(h,'linestyle','-');
set(h,'Color',color);

ylim([0 0.3]);

axis off

print(f,strcat(label,'.eps'),'-depsc')



end

