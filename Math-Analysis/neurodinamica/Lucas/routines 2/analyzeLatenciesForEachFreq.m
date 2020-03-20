%load variables and create empty matrices
load('latency_maxLikelihood_psth_convolved_binsize1_sigma3','complex_cells','complex_i','latency_allCells_maxLikelihood');
% load('partial_corr_byfit_reviewed','responseVectors');
ncells=numel(complex_cells);
lat_bestcond = nan(ncells,1);
i_bestCondition = nan(ncells,1);
bestSFcond = nan(ncells,1);
bestTFcond = nan(ncells,1);
lat_SFs_atBestTF = nan(ncells,8);
lat_TFs_atBestSF = nan(ncells,10);
SF_columnLabels_octv = -4:3;
SF_columnLabels_lin = 2.^SF_columnLabels_octv;
TF_columnLabels_octv = -4:5;
TF_columnLabels_lin = 2.^TF_columnLabels_octv;

%for each cell get latencies for each Sf at best Tf and vice-versa
for i=1:ncells
    
    filename=strcat('multiple_',char(complex_cells(i)));     
    load(filename,'Sf_octave','Tf_octave');
    
    eval(['responseVector=responseVectors.cell' num2str(complex_i(i)) ';']);
    i_bestCondition(i)=find(responseVector==max(responseVector));
    lat_bestcond(i)=latency_allCells_maxLikelihood(i,i_bestCondition(i));
    
    nSFs=numel(Sf_octave); 
    nTFs=numel(Tf_octave); 
    nconditions=nSFs*nTFs; 
    
    bestSFcond(i)=ceil(i_bestCondition(i)/nTFs);
    
    if rem(i_bestCondition(i),nTFs)==0
        bestTFcond(i)=nTFs;
        SFs_i=((nTFs).*(0:nSFs-1))+bestTFcond(i);
    else bestTFcond(i)=rem(i_bestCondition(i),nTFs);
        SFs_i=((nTFs).*(0:nSFs-1))+bestTFcond(i)+1;
    end
    
    TFs_i=(bestSFcond(i)-1)*nTFs+1:bestSFcond(i)*nTFs;
    
    eval(['SFs_i_atBestTF.complexCell' num2str(i) '=SFs_i;']);
    eval(['TFs_i_atBestSF.complexCell' num2str(i) '=TFs_i;']);
    eval(['SFvalues_octave.complexCell' num2str(i) '=Sf_octave;']);
    eval(['TFvalues_octave.complexCell' num2str(i) '=Tf_octave;']);
    
    sf_column_i=find(SF_columnLabels_octv==Sf_octave(1));
    tf_column_i=find(TF_columnLabels_octv==Tf_octave(1));
    
    lat_SFs_atBestTF(i,sf_column_i:sf_column_i+nSFs-1)=latency_allCells_maxLikelihood(i,SFs_i);
    lat_TFs_atBestSF(i,tf_column_i:tf_column_i+nTFs-1)=latency_allCells_maxLikelihood(i,TFs_i);
    
    clear nSFs nTFs nconditions Sf_octave Tf_octave SFs_i TFs_i responseVector sf_column_i tf_column_i filename
    
end

[r_sf,c_sf]=find(~isnan(lat_SFs_atBestTF));
nonempty_cols_sf=nan(1,size(lat_SFs_atBestTF,2));

for k=1:size(lat_SFs_atBestTF,2);
    nonempty_cols_sf(k)=numel(find(c_sf==k));
end

[r_tf,c_tf]=find(~isnan(lat_TFs_atBestSF));
nonempty_cols_tf=nan(1,size(lat_TFs_atBestSF,2));

for k=1:size(lat_TFs_atBestSF,2);
    nonempty_cols_tf(k)=numel(find(c_tf==k));
end

mean_lat_SFs_atBestTF=nanmean(lat_SFs_atBestTF);
sem_lat_SFs_atBestTF=nanstd(lat_SFs_atBestTF)./sqrt(nonempty_cols_sf);
mean_lat_TFs_atBestSF=nanmean(lat_TFs_atBestSF);
sem_lat_TFs_atBestSF=nanstd(lat_TFs_atBestSF)./sqrt(nonempty_cols_tf);
mean_lat_atBestCondition=nanmean(lat_bestcond);
sem_lat_atBestCondition=nanstd(lat_bestcond)/sqrt(numel(find(~isnan(lat_bestcond))));

clear r_sf c_sf r_tf c_tf k

save /Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/latency/latForEachFreq.mat