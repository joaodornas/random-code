load ('latForEachFreq','complex_cells','SFs_i_atBestTF','SFvalues_octave','lat_SFs_atBestTF','SF_columnLabels_octv');
ncells=numel(complex_cells);
analysis_period=800;
analysis_startTime=201;
window_size=10;
baseline_duration=1000;

pref_SFs_dyn=nan(ncells,analysis_period);
BWs_dyn=nan(ncells,analysis_period);
sigmas_dyn=nan(ncells,analysis_period);
LSFVs_dyn=nan(ncells,analysis_period);
LSFSs_dyn=nan(ncells,analysis_period);
HSFSs_dyn=nan(ncells,analysis_period);
response_amp_dyn=nan(ncells,analysis_period);
fit_pvalue=nan(ncells,analysis_period);
rsquare=nan(ncells,analysis_period);
x_axis=(0.05:0.05:64)';

priebe_1D=fittype('A.*((exp((-((log2(x)-log2(ps))^2)./(2.*(tuningWidth+skew.*(log2(x)-log2(ps))).^2))))-(exp(-1./((skew).^2))))');

for i=[1:4,6:12,14:ncells]
    
    tic
    
    disp(['fitting and calculating BWs & LSFSs of cell #' num2str(i)])
    
    fitted_curves_i=nan(numel(x_axis),analysis_period);

    %get SF conditions and values for cell i
    eval(['SFs_i=SFs_i_atBestTF.complexCell' num2str(i) ';']);
    eval(['SFvalues_octv=SFvalues_octave.complexCell' num2str(i) ';']);

    SFvalues_lin=2.^SFvalues_octv;
    eval(['SFvalues_linScale.complexCell' num2str(i) '=SFvalues_lin;']);
    
    %load and shape psths in matrices where each row is a condition
    load(['psth_' char(complex_cells(i)) '_convolved_sigma3.mat'],'psth');

    SF_psths=zeros(numel(SFs_i),numel(psth.condition1));
    for k=1:numel(SFs_i)
        SF_psths(k,:)=eval(['psth.condition' num2str(SFs_i(k))]);
        mean_baseline=mean(SF_psths(k,1:baseline_duration));
        SF_psths(k,:)=SF_psths(k,:)-mean_baseline;
    end

    clear k psth SFs_i
    
    %latencies for each condition, if nan, then lat=55 (mean latency)
    sf_column_i=find(SF_columnLabels_octv==SFvalues_octv(1));
    lat=lat_SFs_atBestTF(i,sf_column_i:sf_column_i+numel(SFvalues_octv)-1)';
    lat(find(isnan(lat)))=55;
    sf_responses=zeros(size(lat));
    sf_responses_cell=nan(numel(SFvalues_octv),analysis_period);
    sf_response_sem=zeros(size(lat));
    sf_response_sem_cell=nan(numel(SFvalues_octv),analysis_period);

    for k=analysis_startTime:analysis_period+analysis_startTime-1 %slide windows, fit and extract parameters
        
        col_i=k-analysis_startTime+1;
        
        for j=1:numel(lat)
        sf_responses(j)=mean(SF_psths(j,baseline_duration+k+lat(j)-(window_size/2):baseline_duration+k+lat(j)+(window_size/2)));
        sf_response_sem(j)=std(SF_psths(j,baseline_duration+k+lat(j)-(window_size/2):baseline_duration+k+lat(j)+(window_size/2)))/sqrt(window_size);
        end
        
        if max(sf_responses)<0.1
            
            eval(['fittedmodels.cell' num2str(i) '.window' num2str(k) '=nan;']);
            eval(['gofs.cell' num2str(i) '.window' num2str(k) '=nan;']);
            eval(['outputs.cell' num2str(i) '.window' num2str(k) '=nan;']);
            
        else
            
        initial_values=[max(sf_responses) SFvalues_lin(find(sf_responses==max(sf_responses))) 0 1.5];
        lower_limits=[0.1,0.05,-1,0.1]; 
        upper_limits=[max(sf_responses)+20,10,1,8];
        [fittedmodel,gof,output]=fit(SFvalues_lin',sf_responses,priebe_1D,'maxIter',1000,'startpoint',initial_values,'lower',lower_limits,'upper',upper_limits);

        [fit_pvalue(i,col_i)]=fitp(gof,output); %evaluate goodness of fit with f stats
        rsquare(i,col_i)=gof.rsquare;

        if rsquare(i,col_i)>0.7 %accept only fits with r2>0.7
            
            pref_SFs_dyn(i,col_i)=fittedmodel.ps;
            sigmas_dyn(i,col_i)=fittedmodel.tuningWidth;
            response_amp_dyn(i,col_i)=fittedmodel.A;
            fitted_curves_i(:,col_i)=feval(fittedmodel,x_axis);
%             BWs_dyn(i,col_i)=fullWidthAtHalfHeight(x_axis,fitted_curves_i(:,col_i),'octave');
            
            sf=SFvalues_lin(1); A=fittedmodel.A; tuningWidth=fittedmodel.tuningWidth; bestSF=fittedmodel.ps; skew=fittedmodel.skew;
            lowest_SF_resp=A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
            sf=bestSF;
            best_SF_resp=A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
            sf=SFvalues_lin(numel(SFvalues_lin));
            highest_SF_resp=A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
            LSFSs_dyn(i,col_i)=lowest_SF_resp/best_SF_resp;
            HSFSs_dyn(i,col_i)=highest_SF_resp/best_SF_resp;
            clear sf lowest_SF_resp best_SF_resp A tuningWidth bestSF skew

        end

        sf_responses_cell(:,col_i)=sf_responses;
        sf_response_sem_cell(:,col_i)=sf_response_sem;
        eval(['fittedmodels.cell' num2str(i) '.window' num2str(k) '=fittedmodel;']);
        eval(['gofs.cell' num2str(i) '.window' num2str(k) '=gof;']);
        eval(['outputs.cell' num2str(i) '.window' num2str(k) '=output;']);

        end
               
    end

    eval(['sf_responses_dyn.cell' num2str(i) '=sf_responses_cell;']);
    eval(['fitted_curves.cell' num2str(i) '=fitted_curves_i;']);
    clear k j SF_psths SFvalues_lin initial_values sf_responses sf_response_sem SFvalues_octv fittedmodel gof output sf_responses_cell fitted_curves_i mean_baseline upper_limits lower_limits lat sf_column_i
    
    disp(['elapsed time for cell ' num2str(i) ' = ' num2str(toc) ' seconds'])
end

% non_empty_cols_lsfs=nan(1,analysis_period);
% non_empty_cols_hsfs=nan(1,analysis_period);
% non_empty_cols_bw=nan(1,analysis_period);
% non_empty_cols_prefSF=nan(1,analysis_period);
% non_empty_cols_sigma=nan(1,analysis_period);
% non_empty_cols_amp=nan(1,analysis_period);
% 
% for k=1:analysis_period
%     non_empty_cols_lsfs(k)=numel(find(~isnan(LSFSs_dyn(:,k))));
%     non_empty_cols_hsfs(k)=numel(find(~isnan(HSFSs_dyn(:,k))));
%     non_empty_cols_bw(k)=numel(find(~isnan(BWs_dyn(:,k))));
%     non_empty_cols_prefSF(k)=numel(find(~isnan(pref_SFs_dyn(:,k))));
%     non_empty_cols_sigma(k)=numel(find(~isnan(sigmas_dyn(:,k))));
%     non_empty_cols_amp(k)=numel(find(~isnan(response_amp_dyn(:,k))));
% end
% 
% mean_lsfs_dyn=nanmean(LSFSs_dyn);
% sem_lsfs_dyn=nanstd(LSFSs_dyn)./sqrt(non_empty_cols_lsfs);
% mean_hsfs_dyn=nanmean(HSFSs_dyn);
% sem_hsfs_dyn=nanstd(HSFSs_dyn)./sqrt(non_empty_cols_hsfs);
% mean_bw_dyn=nanmean(BWs_dyn);
% sem_bw_dyn=nanstd(BWs_dyn)./sqrt(non_empty_cols_bw);
% mean_prefSF_dyn=nanmean(pref_SFs_dyn);
% sem_preSF_dyn=nanstd(pref_SFs_dyn)./sqrt(non_empty_cols_prefSF);
% mean_sigma_dyn=nanmean(sigmas_dyn);
% sem_sigma_dyn=nanstd(sigmas_dyn)./sqrt(non_empty_cols_sigma);
% mean_amp_dyn=nanmean(response_amp_dyn);
% sem_amp_dyn=nanstd(response_amp_dyn)./sqrt(non_empty_cols_amp);

save(['/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/SF_tuningDynamics_windowSize' num2str(window_size) '_800.mat']);

% for i=[1:4,6:12,14:ncells]
%     tic
%     SFvalues_lin=eval(['SFvalues_linScale.complexCell' num2str(i)]);
%     if i==1
%         for k=68:analysis_period
%             disp(['computing LSFVs of cell #' num2str(i) ', time window ' num2str(k) '...'])
%             if ~isnan(pref_SFs_dyn(i,k))
%                 fittedmodel=eval(['fittedmodels.cell' num2str(i) '.window' num2str(k)]);
%                 try
%                     [LSFVs_dyn(i,k)]=calculate_lsfv(fittedmodel,SFvalues_lin(1),'off');
%                 catch
%                     LSFVs_dyn(i,k)=nan;
%                 end
%             end
%         end
%     else
%         for k=1:analysis_period
%             disp(['computing LSFVs of cell #' num2str(i) ', time window ' num2str(k) '...'])
%             if ~isnan(pref_SFs_dyn(i,k))
%                 fittedmodel=eval(['fittedmodels.cell' num2str(i) '.window' num2str(k)]);
%                 try
%                     [LSFVs_dyn(i,k)]=calculate_lsfv(fittedmodel,SFvalues_lin(1),'off');
%                 catch
%                     LSFVs_dyn(i,k)=nan;
%                 end
%             end
%         end
%     end
%     clear SFvalues_lin fittedmodel
%     save(['/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/SF_tuningDynamics_windowSize' num2str(window_size) '.mat']);
%     disp(['elapsed time for cell ' num2str(i) ' = ' num2str(fixdec(toc/60,1)) ' minutes'])
% end
% 
% non_empty_cols_lsfv=nan(1,analysis_period);
% 
% for k=1:analysis_period
%     non_empty_cols_lsfv(i)=numel(find(~isnan(LSFVs_dyn(:,k))));
% end
% 
% mean_lsfv_dyn=nanmean(LSFVs_dyn);
% sem_lsfv_dyn=nanstd(LSFVs_dyn)./sqrt(non_empty_cols_lsfv);
% 
% save(['/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/SF_tuningDynamics_windowSize' num2str(window_size) '.mat']);