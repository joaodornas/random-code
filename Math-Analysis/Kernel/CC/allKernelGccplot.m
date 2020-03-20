function allKernelGccplot

registro = importdata('cc-44px-2sizes-g-Protocols.txt');

start_time = 0;
end_time = 3000;

nConditions = 16;

for r=1:length(registro)

    disp(registro{r});
    
    protocol = registro{r}(2:end-4);
    
    filepath = strcat('/Users/joaodornas/Documents/_Research/_DATA/Center-Surround/kernel/',protocol,'-kernel-density-function-',int2str(start_time),'-',int2str(end_time));
    
    data = load(strcat(filepath,'.mat'));

    timePoints = data.datakernel.kernel(1).timePoints;
    
    for i=1:nConditions
        
        maxDensity(i) = max(data.datakernel.kernel(i).density);
        
    end
    
    for i=1:nConditions
        
        normDensity(i).data = data.datakernel.kernel(i).density ./ max(maxDensity);
        
    end
    
    H = figure;
    
    for i=1:nConditions
        
        subplot(nConditions,1,i);
        
        plot(timePoints,normDensity(i).data,'r');
        
        hold on;
        
    end
    
    print(H,'-depsc',filepath);
       
end


end

