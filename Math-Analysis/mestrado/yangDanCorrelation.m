
function correlationData = yangDanCorrelation(date,site_index,channel,registro,video_index,start_time,end_time,bin_size,nConditions)

tic

disp('BEGIN');

%%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Spass...');

    Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% LOAD CONDITIONS TRIALS LABELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Conditions Trials Labels...');

    for i=1:nConditions
   
        trials_label(i).label = find(Spass.stimIds == i); 
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set Resolution...');

    spike_times = Spass.spike_times;
    
    spike_times = spike_times./ 32000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% BEGIN CONDITIONS LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Begin Conditions Loop...');

    for i=1:nConditions
       
       disp(strcat('Begin Condition . ',int2str(i))); 
        
       trials_spikes = spike_times(trials_label(i).label,:);
        
       nTrial = size(trials_spikes,1);
       
       for t=1:nTrial
          
           spikes_vector = trials_spikes(t,:);
           
           spikes_vector = spikes_vector.';
       
           spikes_vector = sort(spikes_vector);
       
           spikes_vector = spikes_vector(spikes_vector>0);

           spikes_vector = spikes_vector(spikes_vector>(start_time/1000) & spikes_vector<(end_time/1000));
       
           nBins = (end_time - start_time)/bin_size;
           
           correlationData.Condition(i).Trial(t).spike_train = spikes_vector;

           for k=1:nBins
          
             spikes = length(spikes_vector(spikes_vector>=((k-1)*bin_size/1000 + start_time/1000) & spikes_vector<(k*bin_size/1000 + start_time/1000)));
           
             correlationData.Condition(i).Trial(t).binRate(k) = spikes/(nTrial*bin_size/1000);

           end
           
           correlationData.Condition(i).Trial(t).maxBinRate = max(correlationData.Condition(i).Trial(t).binRate);
           
           correlationData.Condition(i).Trial(t).Events(1) = 0;
           correlationData.Condition(i).Trial(t).NonEvents(1) = 0;
           
           for k=1:nBins
              
               if correlationData.Condition(i).Trial(t).binRate(k) > 0.2*correlationData.Condition(i).Trial(t).maxBinRate
                   
                   correlationData.Condition(i).Trial(t).Events(1) = correlationData.Condition(i).Trial(t).Events(1) + length(spikes_vector(spikes_vector>=((k-1)*bin_size/1000 + start_time/1000) & spikes_vector<(k*bin_size/1000 + start_time/1000)));
                   
               else
                   
                   correlationData.Condition(i).Trial(t).NonEvents(1) = correlationData.Condition(i).Trial(t).NonEvents(1) + length(spikes_vector(spikes_vector>=((k-1)*bin_size/1000 + start_time/1000) & spikes_vector<(k*bin_size/1000 + start_time/1000)));
                   
               end
                   
               
           end
           
           correlationData.Condition(i).Trial(t).binRate(correlationData.Condition(i).Trial(t).binRate<0.2*correlationData.Condition(i).Trial(t).maxBinRate) = 0;
  
       end

        
           
           for k=1:nBins
             
                correlationData.Condition(i).Trial(1).CC(1) = 0;
                correlationData.Condition(i).Trial(2).CC(1) = 0;
                correlationData.Condition(i).Trial(nTrial-1).CC(1) = 0;
                correlationData.Condition(i).Trial(nTrial).CC(1) = 0;
                correlationData.Condition(i).Trial(1).meanBinRateBetweenTrials(k) = ( correlationData.Condition(i).Trial(2).binRate(k) + correlationData.Condition(i).Trial(3).binRate(k) + correlationData.Condition(i).Trial(4).binRate(k) + correlationData.Condition(i).Trial(5).binRate(k) ) / 4 ;
                correlationData.Condition(i).Trial(2).meanBinRateBetweenTrials(k) = ( correlationData.Condition(i).Trial(1).binRate(k) + correlationData.Condition(i).Trial(3).binRate(k) + correlationData.Condition(i).Trial(4).binRate(k) + correlationData.Condition(i).Trial(5).binRate(k) ) / 4 ;
                correlationData.Condition(i).Trial(nTrial-1).meanBinRateBetweenTrials(k) = ( correlationData.Condition(i).Trial(nTrial-4).binRate(k) + correlationData.Condition(i).Trial(nTrial-3).binRate(k) + correlationData.Condition(i).Trial(nTrial-2).binRate(k) + correlationData.Condition(i).Trial(nTrial).binRate(k) ) / 4 ;
                correlationData.Condition(i).Trial(nTrial).meanBinRateBetweenTrials(k) = ( correlationData.Condition(i).Trial(nTrial-4).binRate(k) + correlationData.Condition(i).Trial(nTrial-3).binRate(k) + correlationData.Condition(i).Trial(nTrial-2).binRate(k) + correlationData.Condition(i).Trial(nTrial-1).binRate(k) ) / 4 ;
                
           end 
           

           for t=3:nTrial-2
           
                correlationData.Condition(i).Trial(t).CC(1) = 0;
        
               for k=1:nBins
               
                     correlationData.Condition(i).Trial(t).meanBinRateBetweenTrials(k) = ( correlationData.Condition(i).Trial(t - 2).binRate(k) + correlationData.Condition(i).Trial(t - 1).binRate(k) + correlationData.Condition(i).Trial(t + 1).binRate(k) + correlationData.Condition(i).Trial(t + 2).binRate(k) ) / 4 ;

               end
               
               for k=1:nBins
                   
                    correlationData.Condition(i).Trial(t).CC(1) = correlationData.Condition(i).Trial(t).CC(1) + ( 1 / ( nTrial - 1 ) )*((correlationData.Condition(i).Trial(t).binRate(k) - mean(correlationData.Condition(i).Trial(t).binRate)/std(correlationData.Condition(i).Trial(t).binRate))*((correlationData.Condition(i).Trial(t).meanBinRateBetweenTrials(k) - mean(correlationData.Condition(i).Trial(t).meanBinRateBetweenTrials)/std(correlationData.Condition(i).Trial(t).meanBinRateBetweenTrials))));
              
               end
                
           end

           for k=1:nBins
                   
                correlationData.Condition(i).Trial(1).CC(1) = correlationData.Condition(i).Trial(1).CC(1) + ( 1 / ( nTrial - 1 ) )*((correlationData.Condition(i).Trial(1).binRate(k) - mean(correlationData.Condition(i).Trial(1).binRate)/std(correlationData.Condition(i).Trial(1).binRate))*((correlationData.Condition(i).Trial(1).meanBinRateBetweenTrials(k) - mean(correlationData.Condition(i).Trial(1).meanBinRateBetweenTrials)/std(correlationData.Condition(i).Trial(1).meanBinRateBetweenTrials))));
                correlationData.Condition(i).Trial(2).CC(1) = correlationData.Condition(i).Trial(2).CC(1) + ( 1 / ( nTrial - 1 ) )*((correlationData.Condition(i).Trial(2).binRate(k) - mean(correlationData.Condition(i).Trial(2).binRate)/std(correlationData.Condition(i).Trial(2).binRate))*((correlationData.Condition(i).Trial(2).meanBinRateBetweenTrials(k) - mean(correlationData.Condition(i).Trial(2).meanBinRateBetweenTrials)/std(correlationData.Condition(i).Trial(2).meanBinRateBetweenTrials))));
                correlationData.Condition(i).Trial(nTrial-1).CC(1) = correlationData.Condition(i).Trial(nTrial-1).CC(1) + ( 1 / ( nTrial - 1 ) )*((correlationData.Condition(i).Trial(nTrial-1).binRate(k) - mean(correlationData.Condition(i).Trial(nTrial-1).binRate)/std(correlationData.Condition(i).Trial(nTrial-1).binRate))*((correlationData.Condition(i).Trial(nTrial-1).meanBinRateBetweenTrials(k) - mean(correlationData.Condition(i).Trial(nTrial-1).meanBinRateBetweenTrials)/std(correlationData.Condition(i).Trial(nTrial-1).meanBinRateBetweenTrials))));
                correlationData.Condition(i).Trial(nTrial).CC(1) = correlationData.Condition(i).Trial(nTrial).CC(1) + ( 1 / ( nTrial - 1 ) )*((correlationData.Condition(i).Trial(nTrial).binRate(k) - mean(correlationData.Condition(i).Trial(nTrial).binRate)/std(correlationData.Condition(i).Trial(nTrial).binRate))*((correlationData.Condition(i).Trial(nTrial).meanBinRateBetweenTrials(k) - mean(correlationData.Condition(i).Trial(nTrial).meanBinRateBetweenTrials)/std(correlationData.Condition(i).Trial(nTrial).meanBinRateBetweenTrials))));
              
           end
       
       disp(strcat('End Condition . ',int2str(i)));
       
    end
    
%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save data...');
    
    filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));
    filepath = strcat(filepath,'/','correlation','/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel);
    filepath = strcat(filepath,'-bin_size-',int2str(bin_size),'-correlation');
       
    save(filepath,'correlationData');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for i=1:nConditions
       
    f1 = figure;
    
    nTrialsF1 = size(correlationData.Condition(i).Trial,2);
    
    for t=1:nTrialsF1
        
        CC1(t) = correlationData.Condition(1,i).Trial(1,t).CC(1);
        
    end
    
    plot(1:nTrialsF1,CC1,'r');
    ylabel('CC');
    xlabel('Trial');
    
    print(f1,'-djpeg',strcat(filepath,'-plot-CC-condition-',int2str(i)));
    
%     f2 = figure;
%     
%     nTrialsF2 = size(correlationData.Condition(2).Trial,2);
%     
%     for t=1:nTrialsF2
%         
%         CC2(t) = correlationData.Condition(1,2).Trial(1,t).CC(1);
%         
%     end
%     
%     plot(1:nTrialsF2,CC2,'r');
%     ylabel('CC');
%     xlabel('Trial');
% 
%     print(f2,'-djpeg',strcat(filepath,'-plot-CC-condition-2'));
    
    f3 = figure;
    
    for t=1:nTrialsF1
       
        Events1(t) = correlationData.Condition(i).Trial(t).Events(1);
        
    end
    
    plot(1:nTrialsF1,Events1,'r');
    ylabel('Spikes-Events');
    xlabel('Trial');
    
    print(f3,'-djpeg',strcat(filepath,'-plot-Spikes-Events-condition-',int2str(i)));
    
    f4 = figure;
    
    for t=1:nTrialsF1
       
        NonEvents1(t) = correlationData.Condition(i).Trial(t).NonEvents(1);
        
    end
    
    plot(1:nTrialsF1,NonEvents1,'r');
    ylabel('Spikes-NonEvents');
    xlabel('Trial');
    
    print(f4,'-djpeg',strcat(filepath,'-plot-Spikes-NonEvents-condition-',int2str(i)));
    
   end
   
%     f5 = figure;
%     
%     for t=1:nTrialsF2
%        
%         Events2(t) = correlationData.Condition(2).Trial(t).Events(1);
%         
%     end
%     
%     plot(1:nTrialsF2,Events2,'r');
%     ylabel('Spikes-Events');
%     xlabel('Trial');
%     
%     print(f5,'-djpeg',strcat(filepath,'-plot-Spikes-Events-condition-2'));
%     
%     f6 = figure;
%     
%     for t=1:nTrialsF2
%        
%         NonEvents2(t) = correlationData.Condition(2).Trial(t).NonEvents(1);
%         
%     end
%     
%     plot(1:nTrialsF2,NonEvents2,'r');
%     ylabel('Spikes-NonEvents');
%     xlabel('Trial');
%     
%     print(f6,'-djpeg',strcat(filepath,'-plot-Spikes-NonEvents-condition-2'));
    
    
disp('End Conditions Loop');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

disp('END');
    
toc
    
end