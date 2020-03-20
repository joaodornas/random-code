function rasterAll(registro,nConditions,name,latency_file,WOlatency)


latencies = load(char(latency_file));
latency = latencies.latency_time;


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

    spike_times = Spass.spike_times;
    
    spike_times = spike_times./ 32000;
    
    trials_spikes = spike_times(trials_label(1).labels(:),:);
        
    nTrials_cond_one = size(trials_spikes,1);
    
    clear trials_spikes;
    
    trials_spikes = spike_times(trials_label(2).labels(:),:);
        
    nTrials_cond_two = size(trials_spikes,1);
    
    clear trials_spikes;
    
    trials_spikes = spike_times(trials_label(3).labels(:),:);
        
    nTrials_cond_three = size(trials_spikes,1);
    
    clear trials_spikes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   f = figure;
   for i=1:nConditions
       
       disp(strcat('Begin Condition . ',int2str(i))); 
        
       trials_spikes = spike_times(trials_label(i).labels(:),:);
        
       nTrials = size(trials_spikes,1);
       
       if i == 1
           
           cur_color = 'b';
           
           displacement = 0;
           
       elseif i == 2
           
           cur_color = 'r';
           
           displacement = nTrials_cond_one;
           
       elseif i == 3
           
           cur_color = 'g';
           
           displacement = nTrials_cond_one + nTrials_cond_two;
           
       end
 
       for j=1:nTrials
           
           spike_train = trials_spikes(j,:);
           
           if WOlatency
               
               spike_train = spike_train - latency(i);
               
           end
           
           spike_train = spike_train(spike_train>0);
           
           %spike_train = spike_train(spike_train>=(start_time/1000) & spike_train<=(end_time/1000));
        
           spike_train = sort(spike_train);
           
           for k=1:length(spike_train)
       
                plot([spike_train(k) spike_train(k)],[((j-0.8) + displacement + 0.5) (j + displacement + 0.5)],cur_color);
                
                hold on;
    
           end
    
       end
 
   end
   
   xlabel('Tempo (s)');
   ylabel('{Repeti\c{c}\~oes}','interpreter','latex');
   set(gca,'ylim',[1 (nTrials_cond_one + nTrials_cond_two + nTrials_cond_three)]);
   set(gca,'ydir','rev');
  
   if WOlatency
       
        print(f,'-depsc',strcat(name,'-rasterplot-WOlatency'));
        
   else
       
       print(f,'-depsc',strcat(name,'-rasterplot-WITHlatency'));
       
   end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set(gca,'xlim',[range(1) range(2)]);
% set(gca,'xtick',unique([get(gca,'xtick') range]));


end

