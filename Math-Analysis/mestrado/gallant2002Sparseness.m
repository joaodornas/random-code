
function gallant2002Sparseness(date,site_index,channel,registro,video_index,start_time,end_time,bin_size,nConditions)

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
        
       nTrials = size(trials_spikes,1);
       
       spikes_vector = reshape(trials_spikes.',[],1);

       spikes_vector = spikes_vector.';
       
       spikes_vector = sort(spikes_vector);
       
       spikes_vector = spikes_vector(spikes_vector>0);

       spikes_vector = spikes_vector(spikes_vector>(start_time/1000) & spikes_vector<(end_time/1000));
       
       nBins = (end_time - start_time)/bin_size;
       
       for k=1:nBins
          
           spikes = length(spikes_vector(spikes_vector>=((k-1)*bin_size/1000 + start_time/1000) & spikes_vector<(k*bin_size/1000 + start_time/1000)));
           
           sparseData(i).binRate(k) = spikes/(nTrials*bin_size/1000);           
           sparseData(i).binRateSq(k) = ( sparseData(i).binRate(k) )^2;  
                     
       end
  
       disp('...Calculate Sparseness');
       
       sparseData(i).gallant2002SelectivityIndex = ( 1 - ( ( sum(sparseData(i).binRate) )^2 / ( nBins * sum(sparseData(i).binRateSq) ) ) ) / ( 1 - ( 1/nBins ) ) ;
       sparseData(i).gallant2002Sparseness = ( 1 - ( mean(sparseData(i).binRate)^2 / ( mean(sparseData(i).binRate)^2 + std(sparseData(i).binRate)^2 ) ) / ( 1 - ( 1 / nBins ) ) ) ;
       sparseData(i).kurtosis = kurtosis(sparseData(i).binRate) - 3;
       
       sparseData(i).trials_spikes = trials_spikes;
       sparseData(i).nTrials = nTrials;
       
       sparseData(i).SparsenessInTime = 0;
       
       for s=1:size(sparseData(i).binRate,2)
       
            sparseData(i).SparsenessInTime = [sparseData(i).SparsenessInTime ( 1 - ( mean(sparseData(i).binRate(1:s))^2 / ( mean(sparseData(i).binRate(1:s))^2 + std(sparseData(i).binRate(1:s))^2 ) ) / ( 1 - ( 1 / nBins ) ) )];
       
       end
           
        disp(strcat('End Condition . ',int2str(i)));

    end  
    
disp('End Conditions Loop');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save data...');
    
    filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));
    filepath = strcat(filepath,'/','sparseness','/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel);
    filepath = strcat(filepath,'-bin_size-',int2str(bin_size),'-sparseness');
       
    save(filepath,'sparseData');
    
    for i=1:nConditions
    
        f(i) = figure;

        plot(1:size(sparseData(i).SparsenessInTime,2),sparseData(i).SparsenessInTime)
        ylabel(strcat('Sparseness - Condition ',int2str(i)));
        xlabel('Time');

        print(f(i),'-djpeg',strcat(filepath,'-SparsenessInTime-condition',int2str(i)));

%         f2 = figure;
% 
%         plot(1:size(sparseData(2).SparsenessInTime,2),sparseData(2).SparsenessInTime)
%         ylabel('Sparseness - Condition 2');
%         xlabel('Time');
% 
%         print(f2,'-djpeg',strcat(filepath,'-SparsenessInTime-condition2'));
    
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');
    
toc
    
end