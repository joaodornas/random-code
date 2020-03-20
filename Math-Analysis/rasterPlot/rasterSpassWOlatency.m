function f = rasterSpassWOlatency(spike_times,stimIds,condition,start_time,end_time,latency_start_time,cur_color,inverter)

latency_start_time = 500;

latency_end_time = latency_start_time + 300;

p = 10000;

bin_size = 1;

labels = find(stimIds == condition); 

spike_times = spike_times./ 32000;

trials = spike_times(labels(:),:);

alltrials = trials;

alltrials = reshape(alltrials.',[],1);

alltrials = alltrials.';

latency_time = latency(alltrials,latency_start_time,latency_end_time,p); 

if condition == 3

    latency_time = 0;

end

clear alltrials;

disp(strcat('Latencia:',num2str(latency_time)));

   f = figure;
        
   nTrials = length(labels);

   for j=1:nTrials

       spike_train = trials(j,:);
       
       a_e = spike_train(spike_train>=(0/1000) & spike_train<(latency_start_time/1000));

       after_ae = spike_train(spike_train>=(latency_start_time/1000) & spike_train<(end_time/1000)) - latency_time;

       spike_train = [];

       spike_train = [a_e, after_ae];

       spike_train = spike_train(spike_train>0);

       spike_train = sort(spike_train);
       
       if (condition == 2) && (inverter == 1)

                    nBins = (end_time - start_time)/bin_size;

                    %disp(strcat('...invert trial .',int2str(k)));

                    spikes_inverted = invertSpikeTrain(spike_train,start_time,nBins,bin_size);
                    spikes = spikes_inverted.spikes;
                    
                    spike_train = [];
                    spike_train = spikes;

       end

       for k=1:length(spike_train)

            plot([spike_train(k) spike_train(k)],[(j-0.8) j],cur_color);

            hold on;

       end
 
   end
   
   xlabel('Tempo (s)');
   ylabel('{Repeti\c{c}\~oes}','interpreter','latex');
   set(gca,'ylim',[1 nTrials]);
   set(gca,'ydir','rev');
   
   if (condition == 2) && (inverter == 1)
       
        print(f,'-depsc',strcat('condicao-',int2str(condition),'-invertido'));
        
   else
       
       print(f,'-depsc',strcat('condicao-',int2str(condition)));
       
   end
       
       

end

