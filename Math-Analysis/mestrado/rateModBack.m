function rateModBack(date,site_index,channel,registro,video_index,start_time,end_time,bin_size,bin_size_resolution,nConditions)

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
       
       trials_spikes = spike_times(trials_label(i).label,:); 
       allTrials(i) = size(trials_spikes,1);
        
    end

    minTrials = min(allTrials);
    
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
       
       if (i == 2)

            disp('...invert spike train .');
            
            nBins_new = (end_time - start_time)/bin_size_resolution;
            spikes_inverted = invertSpikeTrain(spikes_vector,start_time,nBins_new,bin_size_resolution);
            spikes_vector_inverted = [];
            spikes_vector_inverted = spikes_inverted.spikes;
            
            for k=1:nBins
          
                spikes_inverted = length(spikes_vector_inverted(spikes_vector_inverted>=((k-1)*bin_size/1000 + start_time/1000) & spikes_vector_inverted<(k*bin_size/1000 + start_time/1000)));
           
                rateModBack.conditions(4).binRate(k) = spikes_inverted/(nTrials*bin_size/1000);

            end
       
            rateModBack.conditions(4).meanFR = mean(rateModBack.conditions(4).binRate);
            
            
           for t=1:minTrials

              spikesTrial = trials_spikes(t,:);

              spikesTrial = sort(spikesTrial);

              spikesTrial = spikesTrial(spikesTrial>0);

              spikesTrial = spikesTrial(spikesTrial>(start_time/1000) & spikesTrial<(end_time/1000));

              nBins = (end_time - start_time)/bin_size;
              
              nBins_new = (end_time - start_time)/bin_size_resolution;
              spikes_inverted = invertSpikeTrain(spikesTrial,start_time,nBins_new,bin_size_resolution);
              spikes_vector_inverted = [];
              spikes_vector_inverted = spikes_inverted.spikes;
            

              for k=1:nBins

                    spikes = length(spikes_vector_inverted(spikes_vector_inverted>=((k-1)*bin_size/1000 + start_time/1000) & spikes_vector_inverted<(k*bin_size/1000 + start_time/1000)));

                    rateModBack.conditions(4).trial(t).binRate(k) = spikes/(bin_size/1000);

              end

           end       

       end
       
       for k=1:nBins
          
           spikes = length(spikes_vector(spikes_vector>=((k-1)*bin_size/1000 + start_time/1000) & spikes_vector<(k*bin_size/1000 + start_time/1000)));
           
           rateModBack.conditions(i).binRate(k) = spikes/(nTrials*bin_size/1000);

       end
       
       rateModBack.conditions(i).meanFR = mean(rateModBack.conditions(i).binRate);
       
       for t=1:minTrials
           
          spikesTrial = trials_spikes(t,:);
          
          spikesTrial = sort(spikesTrial);
       
          spikesTrial = spikesTrial(spikesTrial>0);

          spikesTrial = spikesTrial(spikesTrial>(start_time/1000) & spikesTrial<(end_time/1000));
          
          nBins = (end_time - start_time)/bin_size;
          
          for k=1:nBins
          
                spikes = length(spikesTrial(spikesTrial>=((k-1)*bin_size/1000 + start_time/1000) & spikesTrial<(k*bin_size/1000 + start_time/1000)));
           
                rateModBack.conditions(i).trial(t).binRate(k) = spikes/(bin_size/1000);

          end
           
       end       
       
    end
    
    for k=1:nBins
       
        rateModulationDiretoReverso.binRate(k) = rateModBack.conditions(1).binRate(k) / rateModBack.conditions(2).binRate(k);
        rateModulationDiretoInvertido.binRate(k) = rateModBack.conditions(1).binRate(k) / rateModBack.conditions(4).binRate(k);
        rateModulationDiretoInvertido.Contrast(k) = ( rateModBack.conditions(1).binRate(k) - rateModBack.conditions(4).binRate(k) ) / ( rateModBack.conditions(1).binRate(k) + rateModBack.conditions(4).binRate(k) ) ; 
        
    end
    

    rateModBack.rateModulationDiretoReverso = rateModulationDiretoReverso;
    rateModBack.rateModulationDiretoInvertido = rateModulationDiretoInvertido;
    
    filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));
    %filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));
    
    mkdir(filepath,'rateModulation');
    
    filepath = strcat(filepath,'/','rateModulation','/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel);
    filepath = strcat(filepath,'-rate-modulation');
       
    save(filepath,'rateModBack');
    
    
    
    for k=1:nBins
       
        for t=1:minTrials
            
           ratesDireto(t) = rateModBack.conditions(1).trial(t).binRate(k);
           ratesReverso(t) = rateModBack.conditions(2).trial(t).binRate(k);
           ratesInvertido(t) = rateModBack.conditions(4).trial(t).binRate(k);
            
        end
        
        hlillieDireto(k) = lillietest(ratesDireto);
        hlillieReverso(k) = lillietest(ratesReverso);
        hlillieInvertido(k) = lillietest(ratesInvertido);
        
        [pDiretoReverso(k) hDiretoReverso(k)] = signrank(ratesDireto,ratesReverso);
        [pDiretoInvertido(k) hDiretoInvertido(k)] = signrank(ratesDireto,ratesInvertido);
        
    end
    
    
    
    rateModBack.rateModulationDiretoReverso.binRate(isnan(rateModBack.rateModulationDiretoReverso.binRate)) = 0;
    rateModBack.rateModulationDiretoReverso.binRate(isinf(rateModBack.rateModulationDiretoReverso.binRate)) = 0;
    
    
    w(1) = figure;
    plot(0:10/(nBins-1):10,rateModBack.rateModulationDiretoReverso.binRate);
    ylabel('Taxa Relativa','FontSize',24);
    xlabel('Tempo (s)','FontSize',24);
    xlim([0 10]);
    
    hold on;
    
    for i=1:size(rateModBack.rateModulationDiretoReverso.binRate,2)
           
          if ((hlillieDireto(i) == 1) && (hlillieReverso(i) == 1) && (hDiretoReverso(i) == 1) && (pDiretoReverso(i) < 0.05))
               
              if rateModBack.rateModulationDiretoReverso.binRate(i) < 1;

                left = (i/size(rateModBack.rateModulationDiretoReverso.binRate,2))*10;
                right = (i/size(rateModBack.rateModulationDiretoReverso.binRate,2))*10 + (0.5/size(rateModBack.rateModulationDiretoReverso.binRate,2))*10;
                top = rateModBack.rateModulationDiretoReverso.binRate(i);
                bottom = 0;
                x = [left left right right];
                y = [bottom top top bottom];
                  
                  fill(x,y, 'b','EdgeColor','none');
                  
              else
                  
                left = (i/size(rateModBack.rateModulationDiretoReverso.binRate,2))*10;
                right = (i/size(rateModBack.rateModulationDiretoReverso.binRate,2))*10 + (0.5/size(rateModBack.rateModulationDiretoReverso.binRate,2))*10;
                top = rateModBack.rateModulationDiretoReverso.binRate(i);
                bottom = 1;
                x = [left left right right];
                y = [bottom top top bottom];

                  fill(x,y, 'g','EdgeColor','none');
                  
              end
              
          end
          
     end     
    
    p1 = [1 1];
    p2 = [0 size(rateModBack.rateModulationDiretoReverso.binRate,2)];
    plot([p2(1) p2(2)],[p1(1) p1(1)],'Color','r','LineWidth',2);
    
    print(w(1),'-depsc',strcat(filepath,'-Direto-Reverso'));

    rateModBack.rateModulationDiretoInvertido.binRate(isnan(rateModBack.rateModulationDiretoInvertido.binRate)) = 0;
    rateModBack.rateModulationDiretoInvertido.binRate(isinf(rateModBack.rateModulationDiretoInvertido.binRate)) = 0;
    
        
    w(2) = figure;
    plot(0:10/(nBins-1):10,rateModBack.rateModulationDiretoInvertido.binRate);
    ylabel('Taxa Relativa','FontSize',24);
    xlabel('Tempo (s)','FontSize',24);
    xlim([0 10]);

    hold on;
    
         for i=1:size(rateModBack.rateModulationDiretoInvertido.binRate,2)
            
          if ((hlillieDireto(i) == 1) && (hlillieInvertido(i) == 1) && (hDiretoInvertido(i) == 1) && (pDiretoInvertido(i) < 0.05))
        
              if rateModBack.rateModulationDiretoInvertido.binRate(i) < 1;

                left = i;
                right = i + 0.5;
                top = rateModBack.rateModulationDiretoInvertido.binRate(i);
                bottom = 0;
                x = [left left right right];
                y = [bottom top top bottom];
                  
                  fill(x,y, 'b','EdgeColor','none');
                  
              else
                  
                left = (i/size(rateModBack.rateModulationDiretoInvertido.binRate,2))*10;
                right = (i/size(rateModBack.rateModulationDiretoInvertido.binRate,2))*10 + (0.5/size(rateModBack.rateModulationDiretoInvertido.binRate,2))*10;
                top = rateModBack.rateModulationDiretoInvertido.binRate(i);
                bottom = 1;
                x = [left left right right];
                y = [bottom top top bottom];

                  fill(x,y, 'g','EdgeColor','none');
                  
              end
              
          end
          
         end     

    p1 = [1 1];
    p2 = [0 size(rateModBack.rateModulationDiretoInvertido.binRate,2)];
    plot([p2(1) p2(2)],[p1(1) p1(1)],'Color','r','LineWidth',2);

    print(w(2),'-depsc',strcat(filepath,'-Direto-Invertido'));
    
    rateModBack.rateModulationDiretoInvertido.Contrast(isnan(rateModBack.rateModulationDiretoInvertido.Contrast)) = 0;
    rateModBack.rateModulationDiretoInvertido.Contrast(isinf(rateModBack.rateModulationDiretoInvertido.Contrast)) = 0;
    
    w(3) = figure;
    plot(0:10/(nBins-1):10,rateModBack.rateModulationDiretoInvertido.Contrast);
    xlabel('Tempo (s)','FontSize',24);
    ylabel('Contrast','FontSize',24);
    xlim([0 10]);
    print(w(3),'-depsc',strcat(filepath,'-Direto-Invertido-Contrast'));
    
end

