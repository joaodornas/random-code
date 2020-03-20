function rasterWOlatencyGcc(registro,nConditions,start_time,end_time,same_color)


% latency_start_time = 500;
% 
% latency_end_time = latency_start_time + 300;
% 
% p = 10000;

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

       clear alltrials;
       
       if (i == 1) || (i == 4) || (i == 7) || (i == 10) || (i == 13) || (i == 16)
           
           cur_color = 'b';
           
           if i == 1
               
                displacement = 0;
                
           else
               
               displacement = sum(nTrials(1:i-1));
               
           end
           
       elseif (i == 2) || (i == 5) || (i == 8) || (i == 11) || (i == 14)
           
           cur_color = 'r';
           
           displacement = sum(nTrials(1:i-1));
           
       elseif (i == 3) || (i == 6) || (i == 9) || (i == 12) || (i == 15)
           
           cur_color = 'g';
           
           displacement = sum(nTrials(1:i-1));
           
       end
       
       if same_color == 1, cur_color = 'r'; end
 
       for j=1:nTrials(i)
           
           spike_train = trials_spikes(j,:);
           
%            a_e = spike_train(spike_train>=(0/1000) & spike_train<(latency_start_time/1000));
% 
%            after_ae = spike_train(spike_train>=(latency_start_time/1000) & spike_train<(end_time/1000)) - latency_time;
% 
%            spike_train = [];
% 
%            spike_train = [a_e, after_ae];
                          
           spike_train = spike_train(spike_train>0);
           
           spike_train = spike_train(spike_train>=(start_time/1000) & spike_train<=(end_time/1000));
        
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
   
   if same_color == 1
       
       name = strcat('/Users/joaodornas/Documents/_Research/_DATA/Center-Surround/rasterplot/rasterplotWOlatency/',registro(1:end-4),'-rasterplot-WO-latency-same-color.eps');
       
   else
       
       name = strcat('/Users/joaodornas/Documents/_Research/_DATA/Center-Surround/rasterplot/rasterplotWOlatency/',registro(1:end-4),'-rasterplot-WO-latency.eps');
       
   end
      
   print(f,'-depsc',name);
      
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set(gca,'xlim',[range(1) range(2)]);
% set(gca,'xtick',unique([get(gca,'xtick') range]));


end

