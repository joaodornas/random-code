load('population_data','dtc_cells');
n=size(dtc_cells,1);
oriCV_sum_vector=zeros(n,1);
oriCV_int_vector=zeros(n,1);

for i=1:n
    
    %load variables from DTC fit matlab workspace
    dtc_filename=strcat('DTC_',char(dtc_cells(i)));
    load(dtc_filename,'goodness1','output1','fittedmodel1','meanResponse','direction');
    
    %calculate fit p-value
    [p_value_dtc]=fitp(goodness1,output1);
    
    %calculate oriCV with both methods if fitp<0.05 and only sum otherwise
    if p_value_dtc<0.05
        [oriCV_sum,oriCV_int]=oriCV(meanResponse,direction,fittedmodel1,'both');
        oriCV_sum_vector(i)=oriCV_sum;
        oriCV_int_vector(i)=oriCV_int;
    else 
        [oriCV_sum]=oriCV(meanResponse,direction,fittedmodel1,'sum');
        oriCV_sum_vector(i)=oriCV_sum;
        oriCV_int_vector(i)=NaN;
    end
    
end

save('/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/oriCV.mat')