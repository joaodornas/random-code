function allKernelCC

registro = importdata('cc-44px-2sizes-v-Protocols.txt');

for r=1:length(registro)

    disp(registro{r});

    start_time = 0;
    stimulus_start_time = 500;
    stimulus_end_time = 9500;
    end_time = 10000;
    
    calcKernelCC(registro{r},start_time,end_time);
      
end

end

