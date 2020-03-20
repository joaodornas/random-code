function MyContourCorrelation_v4

[min_x,max_x,min_y,max_y] = getMinMax;

plotAll(min_x,max_x,min_y,max_y);

end

function [min_x,max_x,min_y,max_y] = getMinMax

filename = {'AAL-Track' 'Anatomical-Different-AAL-Track' 'Anatomical-Same-AAL-Track' 'Functional-Different-AAL-Track' 'Functional-Same-AAL-Track' 'Voxels-Same-AAL-Track' 'Voxels-Diff-AAL-Track' 'Voxels-Same-Cluster-Track' ...
    'AAL-Passive' 'Anatomical-Different-AAL-Passive' 'Anatomical-Same-AAL-Passive' 'Functional-Different-AAL-Passive' 'Functional-Same-AAL-Passive' 'Voxels-Same-AAL-Passive' 'Voxels-Diff-AAL-Passive' 'Voxels-Same-Cluster-Passive' ...
    'AAL-Rest' 'Anatomical-Different-AAL-Rest' 'Anatomical-Same-AAL-Rest' 'Functional-Different-AAL-Rest' 'Functional-Same-AAL-Rest' 'Voxels-Same-AAL-Rest' 'Voxels-Diff-AAL-Rest' 'Voxels-Same-Cluster-Rest'};

conditions = {'Track' 'Passive' 'Rest'};

min_y = 0;
max_y = 0.1;

min_x = 0;
max_x = 1;

for iLabel=1:length(filename)
   
        load(strcat('Contour-',filename{iLabel},'-bins.mat'));
        
        if min(x_mid(:)) < min_x; min_x = min(x_mid(:)); end
        if max(x_mid(:)) > max_x & max(x_mid(:)) < Inf; max_x = max(x_mid(:)); end

        % if min(y_mid(:)) < min_y; min_y = min(y_mid(:)); end
        if max(y_mid(:)) > max_y; max_y = max(y_mid(:)); end  
    
end

end

function plotAll(min_x,max_x,min_y,max_y)

label = 'AAL';
plotContour(min_x,max_x,min_y,max_y,label);
 
label = 'Anatomical-Different-AAL';
plotContour(min_x,max_x,min_y,max_y,label);

label = 'Anatomical-Same-AAL';
plotContour(min_x,max_x,min_y,max_y,label);

label = 'Functional-Different-AAL';
plotContour(min_x,max_x,min_y,max_y,label);

label = 'Functional-Same-AAL';
plotContour(min_x,max_x,min_y,max_y,label);

label = 'Voxels-Same-AAL';
plotContour(min_x,max_x,min_y,max_y,label);

label = 'Voxels-Diff-AAL';
plotContour(min_x,max_x,min_y,max_y,label);

label = 'Voxels-Same-Cluster';
plotContour(min_x,max_x,min_y,max_y,label);

end

function plotContour(min_x,max_x,min_y,max_y,resolution_label)

z_criterion = 0.13;

conditions = {'Track' 'Passive' 'Rest'};

for iCond=1:length(conditions)
    
    datalabel = strcat(resolution_label,'-',conditions{iCond});

    load(strcat('Contour-',datalabel,'-bins.mat'));

    f = figure;

    hold on;
    a = area([min_x z_criterion],[max_y max_y]);
    set(a,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none');
    contour(x_mid, y_mid, D_xy);
    colormap('jet');
    %plot([0 1],[0 0.75],'k--');
    plot([0 1],[0 1],'k--');
    %plot([0 (max_x+0.1)],[0 0.1*(max_x+0.1)],'k--');
    axis 'square';
    xlabel('mean');
    ylabel('STD');
    xlim([(min_x) (max_x+0.1)]);
    ylim([0 max_y]);
    title(datalabel);
    my_axis = gca;
    set(my_axis,'YTick',[0 0.1 0.2 0.3 0.4 0.5]);
    set(my_axis,'XTick',[0 0.2 0.6 1.0 1.4 1.8]);

    print(f,'-depsc',strcat('Contour-',datalabel,'.eps'));

end

end



