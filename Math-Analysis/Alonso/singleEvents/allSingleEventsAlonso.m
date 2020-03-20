function allSingleEventsAlonso(OS,start_registro,end_registro,Visible)

registro = importdata('memoryBackwardProtocols.txt');

for r=start_registro:end_registro

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
    
    method = 'per_spike';
    
    singleEventAnalysis(registro{r},start_time,end_time,stimulus_start_time,stimulus_end_time,method,Visible,OS);
    
    method = 'threshold';
    
    singleEventAnalysis(registro{r},start_time,end_time,stimulus_start_time,stimulus_end_time,method,Visible,OS);
    
end

end

