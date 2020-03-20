function plotFCgeneralERFC(mu,label,plotParcelsLabels_idx)

doThePlot(mu,label,plotParcelsLabels_idx);

end

function doThePlot(mu,label,plotParcelsLabels_idx)

    hf = figure('color','w');
    hold on
     
    hm = imagesc(mu);

    clrmp = colormap('jet');
    clrmp(1,:) = [1 1 1];

    min_C = 0;
    max_C = 1;

    % colorbar;
    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);
    
    get(hm)
    hold off;
    axis 'equal';
    axis 'off';
    
if plotParcelsLabels_idx
    
    switch plotParcelsLabels_idx
        
        case 1 %%% plot FC MD758 Clusters labels script
            
            run_FC_MD758_Clusters_labels_script;
            
    end
    
end
       
export_fig(label,'-pdf');

close all
    
end