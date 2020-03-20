function allRasterWOlatency(reg)

start_time = 0;

end_time = 10000;

switch reg
    
    case 1

        registro = importdata('memoryBackwardProtocols.txt');
        
    case 2
        
        registro = importdata('cc-44px-2sizes-v-Protocols.txt');
        
end

nConditions = 3;

for r=1:1

    rasterWOlatency(registro{r},nConditions,start_time,end_time,0);
    
    if r == 1
        
        end_time = 4000; 
    
    else
        
        end_time = 10000;
        
    end
    
    rasterWOlatency(registro{r},nConditions,start_time,end_time,1);
    
    close all;
       
end

end

