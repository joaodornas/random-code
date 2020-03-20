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
    unit=dtc_cells(i);
    filename=char(strcat('DTC_',unit));
    load(filename,'conditionMatrix','responseMatrix');
    if size(conditionMatrix,1)>300
        nconditions=size(conditionMatrix,1)/10;
    else nconditions=size(conditionMatrix,1)/5;
    end
    conditionMatrix=log2(conditionMatrix);
    first_sf=conditionMatrix(1,1);
    last_sf=conditionMatrix(size(conditionMatrix,1),1);
    first_tf=conditionMatrix(1,2);
    last_tf=conditionMatrix(size(conditionMatrix,1),2);
    i_mean=zeros(nconditions,1);
    i_var=zeros(size(i_mean));
    for j=first_sf:last_sf
        sf_i=find(conditionMatrix(:,1)==j);
        m=conditionMatrix(sf_i,:);
        for k=first_tf:last_tf
            tf_i=find(m(:,2)==k);
            s=size(tf_i,1);
            trials=responseMatrix(tf_i(1)+sf_i(1)-1:tf_i(1)+sf_i(1)-2+s);
            i_mean((j-first_sf+1)*(k-first_tf+1))=mean(trials);
            i_var((j-first_sf+1)*(k-first_tf+1))=var(trials);
        end
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