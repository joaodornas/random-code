load('population_data','cells','TF_BW');
band_i_TF=find(~isnan(TF_BW));
band_cells_TF=cells(band_i_TF);
n=size(band_i_TF,1);
LSFV_vector_TF=zeros(n,1);
LSFSfromResp_vector_TF=zeros(n,1);
LSFSfromFit_vector_TF=zeros(n,1);

for i=1:n
    
    %load variables from 1D fit matlab workspace
    tf_filename=strcat('TF_',char(band_cells_TF(i)));
    load(tf_filename,'goodness1','output1','fittedmodel1','mean_TF_Values','meanRates_TemporalFreq');
    lowest_TF=mean_TF_Values(1);
    
    %calculate fit p-value
    [p_value_tf]=fitp(goodness1,output1);
    
    %calculate LSFS from actual response
    LSFSfromResp_vector_TF(i)=meanRates_TemporalFreq(1)/max(meanRates_TemporalFreq);
    
    %calculate LSFS and LSFV from fits using calculate_lsfv function if fit
    %is significant
    if p_value_tf<=0.05
        [lsfv,lsfs]=calculate_lsfv(fittedmodel1,lowest_TF,'off');
        LSFV_vector_TF(i)=lsfv;
        LSFSfromFit_vector_TF(i)=lsfs;
    else
        LSFV_vector_TF(i)=NaN; 
        LSFSfromFit_vector_TF(i)=NaN;
    end
    
    save('/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/lsfv_tf.mat')
    
end