function allKernelForBack

registro = importdata('memoryBackwardProtocols.txt');

for r=1:length(registro)

    disp(registro{r});
    
     if r == 1
        
        start_time = 0;
        stimulus_start_time = 500;
        stimulus_end_time = 3500;
        end_time = 4000;
        
    else
        
        start_time = 0;
        stimulus_start_time = 500;
        stimulus_end_time = 9500;
        end_time = 10000;
        
    end
    
    calcKernel(registro{r},start_time,end_time);
      
end

end

