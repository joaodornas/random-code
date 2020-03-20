% load ('latForEachFreq','complex_cells','i_bestCondition','lat_bestcond');
% ncells=numel(complex_cells);
% duration_gaussFit=zeros(ncells,1);
% duration_peak=zeros(ncells,1);
% analysis_period=250;
% baseline_duration=1000;

for i=57:ncells
    
    if i==5 || i==13 %psth's lacking for these cells

        duration_gaussFit(i)=NaN;
        duration_peak(i)=NaN;

    else
        
        %duration of transient by time from latency to peak * 2
        
        psth_i=eval(['psth_bestCondition_complex.cell' num2str(i) ';']);
        psth_period=psth_i(baseline_duration:(baseline_duration+analysis_period));
        i_peakPsth=find(psth_period==max(psth_period),1,'first');
        
        if i_peakPsth > 150 || i_peakPsth-lat_bestcond(i)<10
            duration_peak(i)=nan;
        else duration_peak(i)=(i_peakPsth-lat_bestcond(i))*2;
        end

        %duration of transient by fitting psth with a gaussians and taking 
        %2 SD's 
        time=1:analysis_period+1;
        fresult=fit(time',psth_period','gauss1');
        
        if fresult.b1>150 || fresult.b1<10 || fresult.c1>100
            duration_gaussFit(i)=nan;
        else duration_gaussFit(i) = fresult.c1;
        end
        
        clear psth psth_i parameters i_peakPsth filename cellname psth_period time fresult

    end
    
end

mean_duration_peak=nanmean(duration_peak)
mean_duration_gaussFit=nanmean(duration_gaussFit)
sem_duration_peak=nanstd(duration_peak)./sqrt(numel(find(~isnan(duration_peak))));
sem_duration_gaussFit=nanstd(duration_gaussFit)./sqrt(numel(find(~isnan(duration_gaussFit))));
max_duration_peak=max(duration_peak);
max_duration_gaussFit=max(duration_gaussFit);

mean_duration_byEye=nanmean(duration_byEye);
sem_duration_byEye=nanstd(duration_byEye)./sqrt(numel(find(~isnan(duration_byEye))));
max_duration_byEye=max(duration_byEye);

save /Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/transientDuration.mat