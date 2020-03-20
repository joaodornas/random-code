function preprocess(registro)

start_time = 500;
end_time = 9500;
bin_size = 30;

%%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Spass...');

    Spass = load(strcat(registro,'.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% LOAD CONDITIONS TRIALS label   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Conditions Trials label...');

    nConditions = 1;

    for i=1:nConditions

        trials_label(i).label = find(Spass.stimIds == i); 

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set Resolution...');

    spike_times = Spass.spike_times;

    spike_times = spike_times./ 32000;

    clear Spass;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trials_spikes = spike_times(trials_label(1).label,:);

nTrials = size(trials_spikes,1);

for i=1:nTrials
   
       trial = trials_spikes(i,:);

       trial = trial(trial>0);

       trial = trial(trial>(start_time/1000) & trial<(end_time/1000));

       nBins = (end_time - start_time)/bin_size;

       for k=1:nBins

            spikes = length(trial(trial>=((k-1)*bin_size/1000 + start_time/1000) & trial<(k*bin_size/1000 + start_time/1000)));

            Trials(i).rate(k) = spikes / (bin_size/1000);

       end
    
end

for i=1:nTrials
    
   dlmwrite('response.txt',Trials(i).rate,'delimiter',' ','-append');
        
   dlmwrite('response.txt',[],'newline','pc','-append');
    
end

end