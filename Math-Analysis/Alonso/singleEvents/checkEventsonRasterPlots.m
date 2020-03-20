function checkEventsonRasterPlots(registro,start_time,end_time,results,folderpath,bar,per_bin_per_spike,p,Visible)

nConditions = 3;

Events = results.Events_Conditions;

name = results.name;

inverter = 1;

rasterWOlatencyWIEvents(registro,Events,nConditions,start_time,end_time,inverter,folderpath,bar,name,per_bin_per_spike,p,Visible);


function rasterWOlatencyWIEvents(registro,Events,nConditions,start_time,end_time,inverter,folderpath,bar,name,per_bin_per_spike,p,Visible)

frame_bin_size = 30;
    
latency_start_time = 500;

latency_end_time = latency_start_time + 300;

stimulus_start_time = 500;

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
   set(gcf,'Visible',Visible);
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
           
            condition = i;
            
            protocolo = registro(2:end-4);
            
            latency_time = getLatencyForBack(protocolo,condition);

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
      
       if (i == 2) && (inverter == 1)
           
           h = 4;
           
       else
           
           h = i;
           
       end
       
       nEvents = max(size(Events(h).eventos,1),size(Events(h).eventos,2)); 
       
       if nEvents > 0
           
           for n=1:nEvents

               
               if strcmp(per_bin_per_spike,'per_bin')
                   
                    start = ( Events(h).eventos(n).start * (frame_bin_size/1000) ) + (stimulus_start_time/1000);
                    end_ = ( Events(h).eventos(n).end * (frame_bin_size/1000) ) + (stimulus_start_time/1000); 
               
                    idx =  Events(h).eventos(n).idx;
                    
               elseif strcmp(per_bin_per_spike,'per_spike')
                   
                   start = ( Events(h).eventos(n).start / p ) ;
                    end_ = ( Events(h).eventos(n).end / p ) ; 
                   
                    idx =  Events(h).eventos(n).idx;
                    
               end
                
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
                
                burst = fill(X,Y,'y');
                alpha(burst,0);
                hold on;

           end
           
       end
 
   end
   
   xlabel('Tempo (s)');
   ylabel('{Repeti\c{c}\~oes}','interpreter','latex');
   set(gca,'ylim',[1 (nTrials(1) + nTrials(2) + nTrials(3))]);
   set(gca,'ydir','rev');
  
   
  if (inverter == 1)
   
      print(f,'-depsc',strcat(folderpath,name,bar,name,'-rasterplot-WO-latency-Events-invertido.eps'));
  
  else
      
      print(f,'-depsc',strcat(folderpath,name,bar,name,'-rasterplot-WO-latency-Events.eps'));
      
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set(gca,'xlim',[range(1) range(2)]);
% set(gca,'xtick',unique([get(gca,'xtick') range]));


end


end

