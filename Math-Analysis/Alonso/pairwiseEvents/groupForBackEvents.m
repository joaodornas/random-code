function groupForBackEvents(OS)

if strcmp(OS,'Win')
    
    folderpath = strcat('Z:\DATA\Forward-Backward\Jerome\Alonso\pairwiseEvents\');
    bar = '\';
    
elseif strcmp(OS,'Mac')
    
    folderpath = strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/Alonso/pairwiseEvents/');
    bar = '/';
    
end

registro = importdata('memoryBackwardProtocols.txt');

for r=1:length(registro)
    
  protocol = registro{r}(2:end-4);
  
  database = load(strcat(protocol,'-singleEventsAutoCorr-PSTH-Spike-autocorrelation.mat'));  
  
  Eventos = database.results.Events_Conditions;
  
  Pairwise_Events = groupEvents(Eventos);
   
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
  
  inverter = 1;
   
  nConditions = 2;
  
  mkdir(folderpath,protocol);
  
  rasterWOlatencyPairWiseEvents(registro{r},Pairwise_Events,nConditions,start_time,end_time,stimulus_start_time,stimulus_end_time,inverter,protocol,bar,folderpath);
  
  save(strcat(folderpath,protocol,bar,protocol,'-pairwiseEvents.mat'),'Pairwise_Events');

end





function Pairwise_Events = groupEvents(Events)

    nEvents_For = max(size(Events(1).eventos,1),size(Events(1).eventos,2));
    
    nEvents_Back = max(size(Events(4).eventos,1),size(Events(4).eventos,2));
    
    Pairwise_Events = [];
    
    p = 1;
    for n=1:nEvents_For
        
        idx = Events(1).eventos(n).idx;
        
        start_for = Events(1).eventos(n).start;
        end_for = Events(1).eventos(n).end;
        limiar = Events(1).eventos(n).limiar;
        averageEventTime = Events(1).eventos(n).averageEventTime;
        EventTimeVariability = Events(1).eventos(n).EventTimeVariability;
        
        overlap = false;
        
        e = 1;
        pair = [];
        for m=1:nEvents_Back
            
            start_back = Events(4).eventos(m).start;
            end_back = Events(4).eventos(m).end;
            
            if (start_for < start_back) && (end_for > start_back), overlap = true; end
            if (start_for < end_back) && (end_for > end_back), overlap = true; end
            
            if (start_back < start_for) && (end_back > start_for), overlap = true; end
            if (start_back < end_for) && (end_back > end_for), overlap = true; end
            
            if (start_for == start_back) && (end_for == end_back), overlap = true; end
            
            if overlap == true
                
                pair(e) = m;
                e = e + 1;
                
            end
            
            overlap = false;
            
        end
       
        if ~isempty(pair)
            
            Pairwise_Events(p).Forward_n = n;
            Pairwise_Events(p).Forward_idx = idx;
            Pairwise_Events(p).Forward_start = start_for;
            Pairwise_Events(p).Forward_end = end_for;
            Pairwise_Events(p).Forward_limiar = limiar;
            Pairwise_Events(p).Forward_averageEventTime = averageEventTime;
            Pairwise_Events(p).Forward_EventTimeVariability = EventTimeVariability;
            Pairwise_Events(p).Backward_pair = pair;
            
            if mod(p,2) == 0
                
                cur_color = 'y';
                
            else
                
                cur_color = 'b';
              
            end
            
            Pairwise_Events(p).Forward_color = cur_color;
            Pairwise_Events(p).Backward_color = cur_color;
            
            for e=1:length(pair)
                
                idx_back = Events(4).eventos(pair(e)).idx;
                start_back = Events(4).eventos(pair(e)).start;
                end_back = Events(4).eventos(pair(e)).end;
                limiar_back = Events(4).eventos(pair(e)).limiar;
                averageEventTime_back = Events(4).eventos(pair(e)).averageEventTime;
                EventTimeVariability_back = Events(4).eventos(pair(e)).EventTimeVariability;
                
                Pairwise_Events(p).Backward(e).n = pair(e);
                Pairwise_Events(p).Backward(e).idx = idx_back;
                Pairwise_Events(p).Backward(e).start = start_back;
                Pairwise_Events(p).Backward(e).end = end_back;
                Pairwise_Events(p).Backward(e).limiar = limiar_back;
                Pairwise_Events(p).Backward(e).averageEventTime = averageEventTime_back;
                Pairwise_Events(p).Backward(e).EventTimeVariability = EventTimeVariability_back;
                
            end
            
            p = p + 1;
                
        end
 
    end

end


function rasterWOlatencyPairWiseEvents(registro,Pairwise_Events,nConditions,start_time,end_time,stimulus_start_time,stimulus_end_time,inverter,protocol,bar,folderpath)

frame_bin_size = 30;
    
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
      
       condition = i;
            
       latency_time = getLatencyForBack(protocol,condition);

       clear alltrials;

       disp(strcat('Latencia:',num2str(latency_time)));
       
       if i == 1
           
           cur_color = 'b';
           
           displacement = 0;
           
       elseif i == 2
           
           cur_color = 'r';
           
           displacement = nTrials(1);
           
       end
 
       for j=1:nTrials(i)
           
           spike_train = trials_spikes(j,:);
           
           a_e = spike_train(spike_train>=(start_time/1000) & spike_train<(latency_start_time/1000));

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
       
       if i == 1
           
           nEvents = max(size(Pairwise_Events,1),size(Pairwise_Events,2)); 
       
           if nEvents > 0

               for n=1:nEvents

                    start = ( Pairwise_Events(n).Forward_start * (frame_bin_size/1000) ) + (stimulus_start_time/1000);
                    end_ = ( Pairwise_Events(n).Forward_end * (frame_bin_size/1000) ) + (stimulus_start_time/1000); 

                    idx =  Pairwise_Events(n).Forward_idx;

                    stepX = 2*(frame_bin_size/1000);
                    stepY = 10;

                    p1 = [displacement (nTrials(i) + displacement)];
                    p2 = [start start];
                    line1 = plot([p2(1) p2(1)],[p1(1) p1(2)]);
                    set(line1,'Color','k','LineWidth',0.5);
                    hold on;

                    mytext = text(start + stepX, p1(1) + stepY, int2str(n));
                    set(mytext,'Color','k','FontSize',7);
                    hold on;

                    p1 = [displacement (nTrials(i) + displacement)];
                    p3 = [end_ end_];
                    line2 = plot([p3(1) p3(1)],[p1(1) p1(2)]);
                    set(line2,'Color','k','LineWidth',0.5);
                    hold on;

                    Y = [displacement (nTrials(i) + displacement) (nTrials(i) + displacement) displacement];
                    X = [start start end_ end_];

                    cur_color = Pairwise_Events(n).Forward_color;
                    
                    if ~isempty(Pairwise_Events(n).Backward_pair)
               
                        burst = fill(X,Y,cur_color);
                        hold on;
                        
                    else
                        
                         burst = fill(X,Y,cur_color);
                         alpha(burst,0);
                         hold on;
                         
                    end


               end

           end
       
       elseif i == 2
           
           nEvents = max(size(Pairwise_Events,1),size(Pairwise_Events,2)); 
       
           if nEvents > 0

               for n=1:nEvents
                   
                    if ~isempty(Pairwise_Events(n).Backward_pair)

                        nE = max(size(Pairwise_Events(n).Backward_pair,1),size(Pairwise_Events(n).Backward_pair,2));
                        
                        for e=1:nE
                            
                            start = ( Pairwise_Events(n).Backward(e).start * (frame_bin_size/1000) ) + (stimulus_start_time/1000);
                            end_ = ( Pairwise_Events(n).Backward(e).end * (frame_bin_size/1000) ) + (stimulus_start_time/1000); 

                            idx =  Pairwise_Events(n).Backward(e).idx;

                            stepX = 2*(frame_bin_size/1000);
                            stepY = 10;

                            p1 = [displacement (nTrials(i) + displacement)];
                            p2 = [start start];
                            line1 = plot([p2(1) p2(1)],[p1(1) p1(2)]);
                            set(line1,'Color','k','LineWidth',0.5);
                            hold on;

                            mytext = text(start + stepX, p1(1) + stepY, int2str(n));
                            set(mytext,'Color','k','FontSize',7);
                            hold on;

                            p1 = [displacement (nTrials(i) + displacement)];
                            p3 = [end_ end_];
                            line2 = plot([p3(1) p3(1)],[p1(1) p1(2)]);
                            set(line2,'Color','k','LineWidth',0.5);
                            hold on;

                            Y = [displacement (nTrials(i) + displacement) (nTrials(i) + displacement) displacement];
                            X = [start start end_ end_];
                           
                            cur_color = Pairwise_Events(n).Backward_color;
                            
                            burst = fill(X,Y,cur_color);
                            %alpha(burst,0);
                            hold on;
                        
                        end
                    
                    end
                    
               end
               
           end
          
       end
 
   end
   
   xlabel('Tempo (s)');
   ylabel('{Repeti\c{c}\~oes}','interpreter','latex');
   set(gca,'ylim',[1 (nTrials(1) + nTrials(2))]);
   set(gca,'ydir','rev');
   
   print(f,'-depsc',strcat(folderpath,protocol,bar,protocol,'-rasterplot-WO-latency-Events-invertido.eps'));
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set(gca,'xlim',[range(1) range(2)]);
% set(gca,'xtick',unique([get(gca,'xtick') range]));


end

end

