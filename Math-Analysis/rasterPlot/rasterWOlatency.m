function rasterWOlatency(registro,nConditions,start_time,end_time,inverter)


latency_start_time = 500;

latency_end_time = latency_start_time + 300;

p = 10000;

%%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Spass...');

    Spass = load(char(registro));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% LOAD CONDITIONS TRIALS LABELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Conditions Trials Labels...');

    for i=1:nConditions
   
        trials_label(i).labels = find(Spass.stimIds == i); 
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set Resolution...');

    spike_times = Spass.spike_times ./ 32000;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   f = figure;
   for i=1:nConditions
       
       disp(strcat('Begin Condition . ',int2str(i))); 
        
       trials_spikes = spike_times(trials_label(i).labels(:),:);
        
       nTrials(i) = size(trials_spikes,1);
       
       alltrials = trials_spikes;

       alltrials = reshape(alltrials.',[],1);

       alltrials = alltrials.'; 

       if i == 3

            latency_time = 0;
            
       else
            
            latency_time = latency(alltrials,latency_start_time,latency_end_time,p);

       end

       clear alltrials;

       disp(strcat('Latencia:',num2str(latency_time)));
       
       if i == 1
           
           cur_color = 'b';
           
           displacement = 0;
           
       elseif i == 2
           
           cur_color = 'r';
           
           displacement = nTrials(1);
           
       elseif i == 3
           
           cur_color = 'g';
           
           displacement = nTrials(1) + nTrials(2);
           
       end
 
       for j=1:nTrials(i)
           
           spike_train = trials_spikes(j,:);
           
           a_e = spike_train(spike_train>=(0/1000) & spike_train<(latency_start_time/1000));

           after_ae = spike_train(spike_train>=(latency_start_time/1000) & spike_train<(end_time/1000)) - latency_time;

           spike_train = [];

           spike_train = [a_e, after_ae];
                          
           spike_train = spike_train(spike_train>0);
           
           spike_train = spike_train(spike_train>=(start_time/1000) & spike_train<=(end_time/1000));
        
           spike_train = sort(spike_train);
           
           
           if (i == 2) && (inverter == 1)

                    bin_size  = 1;
               
                    nBins = (end_time - start_time)/bin_size;

                    %disp(strcat('...invert trial .',int2str(k)));

                    spikes_inverted = invertSpikeTrain(spike_train,start_time,nBins,bin_size);
                    spikes = spikes_inverted.spikes;
                    
                    spike_train = [];
                    spike_train = spikes;

           end
           
           for k=1:length(spike_train)
       
                plot([spike_train(k) spike_train(k)],[((j-0.8) + displacement + 0.5) (j + displacement + 0.5)],cur_color);
                
                hold on;
    
           end
    
       end
 
   end
   
   xlabel('Tempo (s)');
   ylabel('{Repeti\c{c}\~oes}','interpreter','latex');
   set(gca,'ylim',[1 (nTrials(1) + nTrials(2) + nTrials(3))]);
   set(gca,'ydir','rev');
  
  if (inverter == 1)
   
      print(f,'-depsc',strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/rasterplot/rasterplotWOlatency-invertido/',registro(1:end-3),'-rasterplot-WO-latency-invertido.eps'));
  
  else
      
      print(f,'-depsc',strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/rasterplot/rasterplotWOlatency/',registro(1:end-3),'-rasterplot-WO-latency.eps'));
      
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set(gca,'xlim',[range(1) range(2)]);
% set(gca,'xtick',unique([get(gca,'xtick') range]));


end

