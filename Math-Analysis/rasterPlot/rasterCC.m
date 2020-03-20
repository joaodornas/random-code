function rasterCC(registro)


%%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Spass...');

    Spass = load(char(registro));
    
    nConditions = Spass.parameters.nconditions;
    
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
    
    for i=1:nConditions
    
        trials_spikes(i).all_spikes = spike_times(trials_label(i).labels(:),:);
        
        nTrials(i) = min(size(trials_spikes(i).all_spikes,1),size(trials_spikes(i).all_spikes,2));
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   f = figure;
   for i=1:nConditions
       
       disp(strcat('Begin Condition . ',int2str(i))); 
       
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
           
           spike_train = trials_spikes(i).all_spikes(j,:);
           
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
   set(gca,'ylim',[1 sum(nTrials(1:end))]);
   set(gca,'ydir','rev');
  
   filepath = strcat('/Users/joaodornas/Documents/_Research/_DATA/Center-Surround/rasterplot/');
   
   print(f,'-depsc',strcat(filepath,registro(1:end-4),'-rasterplot'));
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set(gca,'xlim',[range(1) range(2)]);
% set(gca,'xtick',unique([get(gca,'xtick') range]));

end

