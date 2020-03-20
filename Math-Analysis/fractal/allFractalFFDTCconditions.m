function allFractalFFDTCconditions

tic

registro = importdata('memoryBackwardProtocols.txt');

for r=1:length(registro)
    
   [start_time, end_time, blk] = DTC_forback(char(registro{r}));
      
   if strcmp(blk,'none')
        
        Spass = NaN;
        
   else
        
        Spass = load(char(blk));
        
        getFractal(Spass,start_time,end_time);
    
   end

end


function getFractal(Spass,start_time,end_time)

    nConditions = 16;

    bins = 1000;

    protocol = char(Spass.cellname);

    for n=1:nConditions

        FF_across_bins = load(char(strcat(protocol,'-FF-across_bins-ForBack-DTC-conditions-',int2str(n),'.mat')));

        fractalFanoDTCconditions(FF_across_bins,'DTC','across_bins',n);

        close all;

        clear FF_acroos_bins;

        FF_across_trials = load(char(strcat(protocol,'-FF-across_trials-ForBack-DTC-conditions-',int2str(n),'.mat')));

        fractalFanoDTCconditions(FF_across_trials,'DTC','across_trials',n);

        close all;

        clear FF_across_trials;


    end    
     
end
        
end

