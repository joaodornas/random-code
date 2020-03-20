load('population_data','cells','SF_BW');
band_i=find(~isnan(SF_BW));
band_cells=cells(band_i);
n=size(band_i,1);
LSFV_vector=zeros(n,1);
LSFSfromResp_vector=zeros(n,1);
LSFSfromFit_vector=zeros(n,1);

for i=1:n
    
    %load variables from 1D fit matlab workspace
    sf_filename=strcat('SF_',char(band_cells(i)));
    load(sf_filename,'goodness1','output1','fittedmodel1','mean_SF_Values','meanRates_SpatialFreq');
    lowest_SF=mean_SF_Values(1);
    
    %calculate fit p-value
    [p_value_sf]=fitp(goodness1,output1);
    
    %calculate LSFS from actual response
    LSFSfromResp_vector(i)=meanRates_SpatialFreq(1)/max(meanRates_SpatialFreq);
    
    %calculate LSFS and LSFV from fits using calculate_lsfv function if fit
    %is significant
    if p_value_sf<=0.05
        [lsfv,lsfs]=calculate_lsfv(fittedmodel1,lowest_SF,'off');
        LSFV_vector(i)=lsfv;
        LSFSfromFit_vector(i)=lsfs;
    else
        LSFV_vector(i)=NaN; 
        LSFSfromFit_vector(i)=NaN;
    end
    
    %save('/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/lsfv.mat')
    
end