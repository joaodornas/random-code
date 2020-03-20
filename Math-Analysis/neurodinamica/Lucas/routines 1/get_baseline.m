load('population_data','cells');
n=size(cells,1);
baseline_SEM=zeros(n,1);

for i=1:n
    sf_filename=strcat('SF_',char(cells(i)));
    load(sf_filename,'twoSEM');
    baseline_SEM(i)=twoSEM;
end