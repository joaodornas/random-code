% load ('latForEachFreq','complex_cells','i_bestCondition','lat_bestcond','psth_bestCondition_complex','complex_i');
% load ('transientDuration','duration_peak');
% ncells=numel(complex_cells);
% i_peakPsth=zeros(ncells,1);
% Rmax_vector=zeros(ncells,1);
% Rmin_vector=zeros(ncells,1);
% tau_vector=zeros(ncells,1);
% baseline_vector=zeros(ncells,1);
% fit_pvalue_vector=zeros(ncells,1);
% analysis_period=400;
% baseline_duration=1000;
% stim_duration=4000;

for i=79:ncells
    
    if i==5 || i==13 %psths lacking for these cells

        i_peakPsth(i)=nan;
        Rmax_vector(i)=nan;
        Rmin_vector(i)=nan;
        tau_vector(i)=nan;
        fit_pvalue_vector(i)=nan;
        baseline_vector(i)=nan;

    else
        
        %find initial peak, with limit of 200 ms
        psth_i=eval(['psth_bestCondition_complex.cell' num2str(i) ';']);
        psth_period1=psth_i(baseline_duration+1:(baseline_duration+analysis_period));
        i_peakPsth(i)=find(psth_period1==max(psth_period1),1,'first');

%         if i_peakPsth(i) > 300
%             i_peakPsth(i)=nan;
%         end
%          
        %fit psth from initial peak to end of stimulation period with
        %exponential decay model, as in M?ller et al. (2001) J Neurosci
        %21(17):6978
        if ~isnan(i_peakPsth(i))

            psth_period2=psth_i(baseline_duration+i_peakPsth(i):baseline_duration+stim_duration);
            time=1:numel(psth_period2);
            [fittedmodel,gof,output]=fit(time',psth_period2','(Rmax-Rmin).*exp(-x./tau)+Rmin','startpoint',[max(psth_period2) min(psth_period2) 50],'upper',max(psth_period2));

            eval(['fittedmodels.cell' num2str(i) '=fittedmodel;']);
            eval(['gof.cell' num2str(i) '=gof;']);
            eval(['output.cell' num2str(i) '=output;']);
            eval(['PSTHs_fromPeak.cell' num2str(i) '=psth_period2;']);
            
            %evaluate goodness of fit with f-statistics with alpha=0.05
            [fit_pvalue]=fitp(gof,output);
            
            fit_pvalue_vector(i)=fit_pvalue;
            
            if fit_pvalue >= 0.05
                Rmax_vector(i)=nan;
                Rmin_vector(i)=nan;
                tau_vector(i)=nan;
            else
                Rmax_vector(i)=fittedmodel.Rmax;
                Rmin_vector(i)=fittedmodel.Rmin;
                tau_vector(i)=fittedmodel.tau;
            end

        end
    
        baseline_vector(i)=mean(psth_i(1:baseline_duration));
    
        fig=figure;
        plot(psth_period2);
        hold on
        plot(fittedmodel);
        saveas(fig,['/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/exp_fits_psth_figs/exp_fit_psth_complex_cell' num2str(i) '.fig']);
        close
        clear psth_i psth_period1 psth_period2 Rmax Rmin tau time gof output fit_pvalue fittedmodel
    
    end
    
end

sust_ratio=(Rmin_vector-baseline_vector)./(Rmax_vector-baseline_vector);
trans_ratio=(Rmax_vector-baseline_vector)./(Rmin_vector-baseline_vector);
transientDuration_combined=duration_peak./2+tau_vector;
mean_sust_ratio=nanmean(sust_ratio);
sem_sust_ratio=nanstd(sust_ratio)./sqrt(numel(find(~isnan(sust_ratio))));
mean_trans_ratio=nanmean(trans_ratio);
sem_trans_ratio=nanstd(trans_ratio)./sqrt(numel(find(~isnan(trans_ratio))));
mean_transientDuration_combined=nanmean(transientDuration_combined);
sem_transientDuration_combined=nanstd(transientDuration_combined)./sqrt(numel(find(~isnan(transientDuration_combined))));

save /Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/exp_fit_psth.mat