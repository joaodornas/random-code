cd '\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\others'

%load workspaces containing names of elligible units

load('population_data','x2D_fitted_cells','cells','direction_class')
a=strmatch('directional',direction_class);
b=strmatch('bidirectional',direction_class);
cells=[cells(a);cells(b)];
%n=size(cells,1);
fitted_cells1=x2D_fitted_cells(1:6);
fitted_cells2=x2D_fitted_cells(8:91);
fitted_cells=[fitted_cells1;fitted_cells2];
n=size(fitted_cells,1);

%get sustained cell workspaces and write mean response and frequency
%matrices
condition_mean=[ ];
condition_var=[ ];

for i=1:n
    unit=fitted_cells(i);
    filename=char(strcat('2Dfits_',unit,'.mat'));
    load(filename,'allTrialsResponse','allTrialsDirection');
    a=isnan(allTrialsResponse);
    b=find(a==0);
    responseMatrix=allTrialsResponse(b)';
    conditionMatrix=allTrialsDirection(b)';
    i_mean=zeros(nconditions,1);
    i_var=zeros(size(i_mean));
    step=deg2rad(22.5);
    last=deg2rad(337.5);
    for j=0:step:last
        dir_i=find(conditionMatrix(:,1)==j);
        trials=responseMatrix(dir_i);
        i_mean(j/step+1)=mean(trials);
        i_var(j/step+1)=var(trials);
    end
    condition_mean=[condition_mean;i_mean];
    condition_var=[condition_var;i_var];
    clear('conditionMatrix','responseMatrix');
end

[r,p]=corr(condition_mean,condition_var);
norm_condition_mean=condition_mean/max(condition_mean);
norm_condition_var=condition_var/max(condition_var);
fig=figure;
title 'meanXvar'
plot(norm_condition_mean,norm_condition_var)
save('meanXvar')
saveas(fig,'meanXvar')        