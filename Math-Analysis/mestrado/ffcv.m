function ffcv(date,site_index,channel,registro,start_time,end_time,video_index,bin_size,nConditions)

tic

disp('BEGIN');

%%% READ TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Read trials...');

    Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));
    
    for i=1:nConditions
   
        trials_labels(i).label = find(Spass.stimIds == i); 
    
    end
    
    spike_times = Spass.spike_times./ 32000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');

toc

    allConditionsCV = 0;
    
    for i=1:nConditions
   
        condicoes(i).trials_spikes = spike_times(trials_labels(i).label,:);

        nTrials(i) = size(condicoes(i).trials_spikes,1);
    
        allTrialsCV = 0;
        
        for k=1:nTrials(i)

            spikes = condicoes(i).trials_spikes(k,:);
            spikes = spikes(spikes>0);
            spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));

            condicoes(i).trials(k).spikes = sort(spikes);
            nSpikes = length(condicoes(i).trials(k).spikes);
            
%             if operational_time == 1
%                 
%                 condicoes(i).trials(k).optBin = sshist(condicoes(i).trials(k).spikes);
%                 [condicoes(i).trials(k).RateDensity condicoes(i).trials(k).RateTime] = ksdensity(condicoes(i).trials(k).spikes,'width',condicoes(i).trials(k).optW);
% 
%             end
            
            condicoes(i).trials(k).distances = [];
            
            for j=1:(nSpikes - 1)
               
                d = condicoes(i).trials(k).spikes(j+1) - condicoes(i).trials(k).spikes(j);
                
                condicoes(i).trials(k).distances = [condicoes(i).trials(k).distances d];                
                
            end
            
            condicoes(i).trials(k).meanDistances = mean(condicoes(i).trials(k).distances);
            condicoes(i).trials(k).stdDistances = std(condicoes(i).trials(k).distances);
            
            condicoes(i).trials(k).coefficientOfVariation = condicoes(i).trials(k).stdDistances/condicoes(i).trials(k).meanDistances;
            
            allTrialsCV = allTrialsCV + condicoes(i).trials(k).coefficientOfVariation;
            
            nBins = (end_time - start_time) / bin_size;
            
            for b=1:nBins
               
                nSpikes = length(spikes(spikes>=((b-1)*bin_size/1000 + start_time/1000) & spikes<(b*bin_size/1000 + start_time/1000)));
                
                condicoes(i).trials(k).binRate(b) = nSpikes;
                
            end
            
            
        end
        
        condicoes(i).meanCVtrials = allTrialsCV / nTrials(i);
        
        allConditionsCV = allConditionsCV + condicoes(i).meanCVtrials;
        
       for b=1:nBins
            
          for t=1:nTrials(i)
             
              rate(t) = condicoes(i).trials(t).binRate(b);
              
          end
          
          condicoes(i).FF(b) = var(rate) / mean(rate) ;
           
       end
        
        
    end
    
    
     meanCVconditions = allConditionsCV / length(condicoes);
     

     %%% SALVA NO DISCO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save data...');

%filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));
filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));

filepath = strcat(filepath,'/','ff-cv','/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel);

filepath = strcat(filepath,'-ff-cv');
    
save(filepath,'condicoes');


for i=1:nConditions

  f(i) = figure;
     
     plot(1:nBins,condicoes(i).FF,'b');
     ylabel('FF');
     xlabel('Time Bins');
     
     print(f(i),'-djpeg',strcat(filepath,'-plot-FF-condition-',int2str(i)));
     
%      f2 = figure;
%      
%      plot(1:nBins,condicoes(2).FF,'b');
%      ylabel('FF');
%      xlabel('Time Bins');
%      
%      print(f2,'-djpeg',strcat(filepath,'-plot-FF-condition-2'));

end

end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%