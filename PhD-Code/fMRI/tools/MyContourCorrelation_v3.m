function MyContourCorrelation_v3

label = 'AAL';
plotContour(label);
 
% label = 'Anatomical-Different-AAL';
% plotContour(label);
% 
% label = 'Anatomical-Same-AAL';
% plotContour(label);

% label = 'Functional-Different-AAL';
% plotContour(label);
% 
% label = 'Functional-Same-AAL';
% plotContour(label);

% label = 'Voxels-Same-AAL';
% plotContour(label);
% 
% label = 'Voxels-Different-AAL';
% plotContour(label);
% 
% label = 'Voxels-Same-Cluster';
% plotContour(label);

end

function plotContour(resolution_label)

conditions = {'Track' 'Passive' 'Rest'};

for iCond=1:length(conditions)
    
    datalabel = strcat(resolution_label,'-',conditions{iCond});

    load(strcat('Contour-',datalabel,'-bins.mat'));

    f = figure;

    hold on;
    contour(x_mid, y_mid, D_xy);
    contourcmap('jet','Colorbar','on');
    plot([0 0.5],[0 0.5],'k-');
    if min(x_mid(:)) < 0; plot([0 -0.5],[0 0.5],'k-'); end
    axis 'square';
    xlabel('mean');
    ylabel('STD');
    if min(x_mid(:)) < 0; xlim([-1 1]); else xlim([0 1]); end
    ylim([0 0.5]);
    title(datalabel);

    print(f,'-depsc',strcat('Contour-',datalabel,'.eps'));

end

end



