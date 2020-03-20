load ('latForEachFreq','complex_cells','i_bestCondition','lat_bestcond','psth_bestCondition_complex','complex_i');

ncells=numel(complex_cells);
transience_index=zeros(ncells,1);
p_value_vector=zeros(ncells,1);
sustORtrans_signrank=num2str(zeros(ncells,1));
analysis_period=50;
baseline_duration=1000;


for i=1:ncells
    
    if isnan(lat_bestcond(i));

        transience_index(i)=nan; 
        p_value_vector(i)=nan;
        sustORtrans_signrank(i)='n';
        
    else
        
        psth_i=eval(['psth_bestCondition_complex.cell' num2str(i) ';']);
        period1=psth_i(baseline_duration+lat_bestcond(i):baseline_duration+lat_bestcond(i)+analysis_period);
        period2=psth_i(baseline_duration+lat_bestcond(i)+analysis_period+100:baseline_duration+lat_bestcond(i)+analysis_period+100+analysis_period);
        
        transience_index(i)=mean(period1)/(mean(period1)+mean(period2));
        [p_value_vector(i)]=signrank(period1,period2);
        
        if p_value_vector(i)>0.05
            sustORtrans_signrank(i)='s';
        else sustORtrans_signrank(i)='t';
        end
        
        clear psth_i period1 period2
    end
end

save /Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/transORsust.mat