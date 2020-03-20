slopes_lat_SFs=nan(92,1);
slopes_lat_TFs=nan(92,1);
x_curve_sf=SF_columnLabels_octv;
x_curve_tf=TF_columnLabels_octv;

for i=[1:4,6:12,14:92]

    try
        curve=lat_SFs_atBestTF(i,:);
        curve=curve(find(~isnan(curve)));
        x=x_curve_sf(find(~isnan(curve)));
        fitResult=fit(x',curve','poly1');
        slopes_lat_SFs(i)=fitResult.p1;
    catch
        slopes_lat_SFs(i)=nan;
    end

    try
        curve=lat_TFs_atBestSF(i,:);
        curve=curve(find(~isnan(curve)));
        x=x_curve_tf(find(~isnan(curve)));
        fitResult=fit(x',curve','poly1');
        slopes_lat_TFs(i)=fitResult.p1;
    catch
        slopes_lat_TFs(i)=nan;
    end

    clear curve fitResult x

end

clear x_curve_sf x_curve_tf

mean_slopes_lat_SFs=nanmean(slopes_lat_SFs)
sem_slopes_lat_SFs=nanstd(slopes_lat_SFs)/sqrt(numel(find(~isnan(slopes_lat_SFs))));
mean_slopes_lat_TFs=nanmean(slopes_lat_TFs)
sem_slopes_lat_TFs=nanstd(slopes_lat_TFs)/sqrt(numel(find(~isnan(slopes_lat_TFs))));

p_signrank_slopes_lat_SFs=signrank(slopes_lat_SFs)
p_signrank_slopes_lat_TFs=signrank(slopes_lat_TFs)