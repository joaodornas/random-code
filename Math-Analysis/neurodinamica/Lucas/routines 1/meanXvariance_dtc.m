cd '\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\others'

%load workspaces containing names of elligible units

load('population_data','dtc_cells')
n=size(dtc_cells,1);

%get sustained cell workspaces and write mean response and frequency
%matrices
condition_mean=[ ];
condition_var=[ ];

for i=1:n
    unit=dtc_cells(i);
    filename=char(strcat('DTC_',unit,'.mat'));
    load(filename,'allTrialsResponse','allTrialsDirection');
    a=isnan(allTrialsResponse);
    b=find(a==0);
    responseMatrix=allTrialsResponse(b)';
    conditionMatrix=allTrialsDirection(b)';
    conditionMatrix=single(rad2deg(conditionMatrix));
    i_mean=zeros(16,1);
    i_var=zeros(size(i_mean));
    for j=0:22.5:337.5
        dir_i=find(conditionMatrix==j);
        trials=responseMatrix(dir_i);
        i_mean(j/22.5+1)=mean(trials);
        i_var(j/22.5+1)=var(trials);
    end
    condition_mean=[condition_mean;i_mean];
    condition_var=[condition_var;i_var];
end

[r,p]=corr(condition_mean,condition_var,'rows','pairwise');
norm_condition_mean=condition_mean/max(condition_mean);
norm_condition_var=condition_var/max(condition_var);
fig=figure;
title 'meanXvar_dtc'
plot(norm_condition_mean,norm_condition_var)
save('meanXvar_dtc')
saveas(fig,'meanXvar_dtc')        