function Events_Conditions = findEventsPerSpikeTime(kernel,Spike_Trains_Conditions,spike_distance,silence_distance,stimulus_start_time,stimulus_end_time,limiar_spikes,p)
    
    for g=1:length(Spike_Trains_Conditions)
       
       eventos = getEvents(kernel(g).density,Spike_Trains_Conditions(g).allTrials,spike_distance,silence_distance,stimulus_start_time,stimulus_end_time,limiar_spikes,p);
       
       Events_Conditions(g).eventos = eventos;
      
    end
    
    
function eventos = getEvents(density,allTrials,spike_distance,silence_distance,stimulus_start_time,stimulus_end_time,limiar_spikes,p)
    
    nTrials = max(size(allTrials,1),size(allTrials,2));
    
    block = [];
    for d=1:nTrials
        
        spikeTrain = allTrials(d).trial;

        block = [block, spikeTrain];
        
    end
    
    TotalSpikes = max(size(block,1),size(block,2));
    
    block = reshape(block.',[],1);
    
    block(block==0) = [];
    
    block = sort(block);
    
    block = block(block>=(stimulus_start_time/1000) & block<=(stimulus_end_time/1000));
    
    block = sort(block);
    
    nSpikes = numel(block);
    
    is_counting = false;
    b = 1;
    number_of_spikes = 0;
    for d=2:(nSpikes-3)
        
        past = block(d-1);
        
        firstSpike = block(d);
        
        nextSpike = block(d+1);
        
        distance_spike_time = nextSpike - firstSpike;
        
        [meanDistanceInPast, stdDistanceInPast] = getMeanSilenceDistanceInPast(firstSpike,allTrials);
        
        [meanDistanceInFuture, stdDistanceInFuture] = getMeanSilenceDistanceInFuture(firstSpike,allTrials);
        
        meanDistanceFromEndTime = abs(firstSpike - (stimulus_end_time/1000));

        if ~is_counting
        
            if (distance_spike_time < (spike_distance/1000)) && (((meanDistanceInPast + stdDistanceInPast) >= (silence_distance/1000) || (meanDistanceInPast - stdDistanceInPast) >= (silence_distance/1000))) 
            
                is_counting = true;
                
                spiketime = block(d);
                
                burst(b).firstSpike = spiketime; 
                
                burst(b).meanDistanceInPast = meanDistanceInPast;
                
                burst(b).stdDistanceInPast = stdDistanceInPast;
                
                burst(b).silence_distance = silence_distance/1000;
                
                number_of_spikes = 1;
      
            end 
            
        else
            
            if (distance_spike_time < (spike_distance/1000))
                
                is_counting = true;
                
                number_of_spikes = number_of_spikes + 1;
                
                if (d+1) == (nSpikes-3)
                    
                        number_of_spikes = number_of_spikes + 1;

                        is_counting = false;

                        spiketime = block(d+1);

                        burst(b).lastSpike = spiketime;
        
                        [meanDistanceInFuture, stdDistanceInFuture] = getMeanSilenceDistanceInFuture(spiketime,allTrials);

                        burst(b).meanDistanceInFuture = meanDistanceInFuture;

                        burst(b).stdDistanceInFuture = stdDistanceInFuture;

                        burst(b).number_of_spikes = number_of_spikes;

                        b = b + 1;

                        number_of_spikes = 0;
                    
                end
                
            else
            
                if (distance_spike_time >= (spike_distance/1000)) && (((meanDistanceInFuture + stdDistanceInFuture) >= (silence_distance/1000) || (meanDistanceInFuture - stdDistanceInFuture) >= (silence_distance/1000)))

                        number_of_spikes = number_of_spikes + 1;

                        is_counting = false;

                        spiketime = block(d);

                        burst(b).lastSpike = spiketime;

                        burst(b).meanDistanceInFuture = meanDistanceInFuture;

                        burst(b).stdDistanceInFuture = stdDistanceInFuture;

                        burst(b).number_of_spikes = number_of_spikes;

                        b = b + 1;

                        number_of_spikes = 0;
 

                elseif (meanDistanceFromEndTime) <= (spike_distance/1000)

                    if number_of_spikes >= limiar_spikes

                        number_of_spikes = number_of_spikes + 1;

                        is_counting = false;

                        spiketime = block(d);

                        burst(b).lastSpike = spiketime;

                        burst(b).meanDistanceInFuture = meanDistanceInFuture;

                        burst(b).stdDistanceInFuture = stdDistanceInFuture;

                        burst(b).number_of_spikes = number_of_spikes;
                        
                        b = b + 1;

                        number_of_spikes = 0;

                    else
                       
                        is_counting = false;

                        number_of_spikes = 0;

                        burst(b) = [];

                    end

                elseif (distance_spike_time >= (spike_distance/1000)) && ((meanDistanceInFuture + stdDistanceInFuture) < (silence_distance/1000) && (meanDistanceInFuture - stdDistanceInFuture) < (silence_distance/1000))

                    if (abs(firstSpike - block(d+2) < (spike_distance/1000)) && abs(firstSpike - block(d+3)) < (spike_distance/1000)) 

                        is_counting = true;

                        number_of_spikes = number_of_spikes + 1;

                    else
                       
                        is_counting = false;

                        number_of_spikes = 0;

                        burst(b) = [];

                    end

                end
            
            end
            
        end
         
    end        

    nBursts = max(size(burst,1),size(burst,2));
    
    if nBursts > 0
        
        if isempty(burst(nBursts).lastSpike)
        
            burst(nBursts).lastSpike = block(nSpikes);
        
        end
        
  
    %     for n=1:nBursts
    %        
    %         firstSpike = burst(n).firstSpike;
    %         lastSpike = burst(n).lastSpike;
    %         
    %         if lastSpike - firstSpike >= spike_distance/1000
    %             
    %             vector(n) = 1;
    %             new_burst(n).firstSpike = firstSpike;
    %             new_burst(n).lastSpike = lastSpike;
    %             
    %         else
    %             
    %             vector(n) = 0;
    %             new_burst(n).firstSpike = firstSpike;
    %             new_burst(n).lastSpike = lastSpike;
    %             
    %         end
    %         
    %     end
    %     
    %     burst = [];
    %     
    %     for n=1:nBursts
    %        
    %         s = 1;
    %         if vector(n) == 1
    %             
    %             burst(s).firstSpike = new_burst(n).firstSpike;
    %             burst(s).lastSpike = new_burst(n).lastSpike;
    %             
    %             s = s + 1;
    %             
    %         end
    %         
    %     end
    %     
    %     clear new_burst;


        for b=1:nBursts

            for d=1:nTrials

                spikes = allTrials(1,d).trial;

                spikes = sort(spikes);

                inferior_limit = burst(b).firstSpike;

                superior_limit = burst(b).lastSpike;

                if numel(spikes) > 0

                    spikes = spikes(spikes>=(inferior_limit) & spikes<=(superior_limit));

                    spikes = sort(spikes);

                    Total_Spikes = numel(spikes);

                    if ~isempty(spikes)

                        firstspike(d) = spikes(1);
                        lastspike(d) = spikes(end);

                    else

                        firstspike(d) = 0;
                        lastspike(d) = 0;

                    end

                end

            end

            firstspike(firstspike==0) = [];
            lastspike(lastspike==0) = [];

            burst(b).averageEventStartTime = mean(firstspike);
            burst(b).EventStartTimeVariability = std(firstspike);

            burst(b).averageEventEndTime = mean(lastspike);
            burst(b).EventEndTimeVariability = std(lastspike);

            burst(b).eventMax = ( TotalSpikes.*max(density( ( (burst(b).firstSpike)*p ) : (burst(b).lastSpike)*p )) ) / nTrials ;

            burst(b).start = burst(b).firstSpike*p;
            burst(b).end = burst(b).lastSpike*p;
            burst(b).idx = find(density==max(density( ( (burst(b).firstSpike)*p ) : (burst(b).lastSpike)*p )));

        end

        s = 1;
        for b=1:nBursts

           if limiar_spikes  <= burst(b).number_of_spikes;

               eventos(s).firstSpike =  burst(b).firstSpike; 

               eventos(s).lastSpike = burst(b).lastSpike;

               eventos(s).meanDistanceInPast = burst(b).meanDistanceInPast;

               eventos(s).stdDistanceInPast = burst(b).stdDistanceInPast;

               eventos(s).meanDistanceInFuture = burst(b).meanDistanceInFuture;

               eventos(s).stdDistanceInFuture = burst(b).stdDistanceInFuture;

               eventos(s).silence_distance = burst(b).silence_distance;

               eventos(s).averageEventStartTime = burst(b).averageEventStartTime;

               eventos(s).EventStartTimeVariability = burst(b).EventStartTimeVariability;

               eventos(s).averageEventEndTime = burst(b).averageEventEndTime;

               eventos(s).EventEndTimeVariability = burst(b).EventEndTimeVariability;

               eventos(s).eventMax = burst(b).eventMax;

               eventos(s).start = burst(b).start;

               eventos(s).end = burst(b).end;

               eventos(s).idx = burst(b).idx;

               eventos(s).durationInterval = burst(b).lastSpike - burst(b).firstSpike;

               eventos(s).number_of_spikes = burst(b).number_of_spikes;

               s = s + 1;

           end

        end
    
    end
    
end
    
    
    

    
end
    
    

% function overlap = eventsOverlap(eventos)
% 
%     
%     nEventos = max(size(eventos,1),size(eventos,2));
%     
%     overlap = false;
%     
%     for n=1:(nEventos - 1)
%                 
%        start_antes = eventos(n).start;
%        end_antes = eventos(n).end;
%        
%        start_depois = eventos(n+1).start;
%        end_depois = eventos(n+1).end;
%        
%        if end_antes >= start_depois
%            
%            overlap = true;
%            
%        end
%                 
%     end
%     
% 
% end

% function grupo = checkEventsIsolated(eventos)
% 
%     nEventos = max(size(eventos,1),size(eventos,2));
% 
%     grupo = [];
%     
%     s = 1;
% 
%     if nEventos > 0
% 
%         for n=2:(nEventos-2)
% 
%             start_primeiro = eventos(n-1).start;
%             end_primeiro = eventos(n-1).end;
% 
%             start_antes = eventos(n).start;
%             start_depois = eventos(n+1).start;
% 
%             end_antes = eventos(n).end;
%             end_depois = eventos(n+1).end;
% 
%             start_segundo = eventos(n+2).start;
%             end_segundo = eventos(n+2).end;
% 
%             if (start_depois - end_antes) == 1
% 
%                 juntos = true;
% 
%             else
% 
%                 juntos = false;
% 
%             end
% 
%             if (start_antes - end_primeiro) > 1
% 
%                 afastado_primeiro = true;
% 
%             else
% 
%                 afastado_primeiro = false;
% 
%             end
% 
%             if (start_segundo - end_depois) > 1
% 
%                 afastado_segundo = true;
% 
%             else
% 
%                 afastado_segundo = false;
% 
%             end
% 
%             if (juntos && afastado_primeiro && afastado_segundo)
% 
%                 grupo(s) = n;
% 
%                 s = s + 1;
% 
%             else
% 
%                 grupo(s) = 0;
% 
%                 s = s + 1;
% 
%             end
% 
%         end
% 
%     else
% 
%         grupo(1) = 0;
% 
%     end
% 
%     
%     if isempty(grupo), grupo(1) = 0; end
%     
% 
% end

% function new_eventos = isolateEventos(rate,eventos,isolados)
% 
% nEventos = max(size(eventos,1),size(eventos,2));
% 
% new_eventos(1) = eventos(1);
% 
% s = 1;
% for n=2:(nEventos-2)
% 
%     idx_anterior = find(isolados==(n+1));
% 
%     if length(idx_anterior) == 0
% 
%         if n == 2, s = s + 1; end
% 
%         new_eventos(s) = eventos(n);
% 
%     else 
% 
%         if n == 2, s = s + 1; end
% 
%         new_eventos(s).start = eventos(n).start;
%         new_eventos(s).end = eventos(n+1).end;
%         
%         start = new_eventos(s).start;
%         end_ = new_eventos(s).end;
%         
%         new_eventos(s).idx = start + find(rate(start:end_)==max(rate(start:end_)));
%         new_eventos(s).limiar = eventos(n).limiar;
%         new_eventos(s).averageEventTime = mean([eventos(n).averageEventTime, eventos(n+1).averageEventTime]);
%         new_eventos(s).EventTimeVariability = mean([eventos(n).EventTimeVariability, eventos(n+1).EventTimeVariability]);
% 
%     end
%     
%     s = s + 1;
% 
% end
% 
% end
