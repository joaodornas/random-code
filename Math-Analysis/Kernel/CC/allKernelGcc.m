function allKernelGcc

registro = importdata('cc-44px-2sizes-g-Protocols.txt');

for r=1:length(registro)

    disp(registro{r});

    start_time = 0;
    stimulus_start_time = 500;
    stimulus_end_time = 2500;
    end_time = 3000;
    
    calcKernelGcc(registro{r},start_time,end_time);
      
end

end

