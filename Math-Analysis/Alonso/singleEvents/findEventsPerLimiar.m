function Events_Conditions = findEventsPerLimiar(PSTH_Conditions,Spike_Trains_Conditions,threshold,stimulus_start_time,frame_bin_size)
    
    for g=1:length(PSTH_Conditions)
       
       limiar = threshold; 
        
       size_rate = max(size(PSTH_Conditions(g).rate,1),size(PSTH_Conditions(g).rate,2));
       
       eventos = getEvents(PSTH_Conditions(g).rate,Spike_Trains_Conditions(g).allTrials,limiar,stimulus_start_time,frame_bin_size,size_rate);
       
       while eventsOverlap(eventos)
           
           limiar = limiar + 2;
           
           eventos = [];
           
           eventos = getEvents(PSTH_Conditions(g).rate,Spike_Trains_Conditions(g).allTrials,limiar,stimulus_start_time,frame_bin_size,size_rate);
           
       end 
       
%       if ~isempty(eventos) 
%           
%             isolados = checkEventsIsolated(eventos);
%     
%             if max(size(isolados,1),size(isolados,2)) > 0
% 
%                 eventos = isolateEventos(PSTH_Conditions(g).rate,eventos,isolados);
% 
%             end
%             
%       end
       
       Events_Conditions(g).eventos = eventos;
      
    end
    
    
function eventos = getEvents(rate,allTrials,limiar,stimulus_start_time,frame_bin_size,size_rate)

     s = 1;
       
       tem_eventos = 'nao';
       
       k = 1;
       
       for h=1:size_rate
          
          if k < size_rate
              
              if rate(k) >= limiar

                   eventos(s).idx = k;

                   eventos(s).start = rateBefore(rate,k);

                   eventos(s).end = rateAfter(rate,k,size_rate);

                   eventos(s).limiar = limiar;

                   fim = eventos(s).end;

                   start = eventos(s).start;

                   nTrials = max(size(allTrials,1),size(allTrials,2));

                   k = fim + 3;
                    
                   for n=1:nTrials

                       spikes = allTrials(1,n).trial;

                       spikes = sort(spikes);

                       if numel(spikes) > 0

                            spikes = spikes(spikes>=((stimulus_start_time + (start-1)*frame_bin_size)/1000) & spikes<=((stimulus_start_time + (fim)*frame_bin_size)/1000));
                            
                            if numel(spikes) > 0
                                
                                spikes = sort(spikes);
                                
                                first_spike = spikes(1);
                                
                            else
                                
                                first_spike = 0;
                                
                            end

                            if first_spike > 0
                                
                                spikeTime(n) = first_spike;
                                
                                
                            else
                                
                                spikeTime(n) = 0;
                                
                            end
                            
                       else

                            spikeTime(n) = 0;

                       end

                   end
                    
                   spikeTime(spikeTime==0) = [];
                   
                   eventos(s).averageEventTime =  mean(spikeTime);

                   eventos(s).EventTimeVariability = std(spikeTime);
                   
                   s = s + 1;

                   tem_eventos = 'sim';

              end

               k = k + 1;
           
          end
          
       end 
    
        if strcmp(tem_eventos,'nao')

            eventos = [];

       end

end

function overlap = eventsOverlap(eventos)

    
    nEventos = max(size(eventos,1),size(eventos,2));
    
    overlap = false;
    
    for n=1:(nEventos - 1)
                
       start_antes = eventos(n).start;
       end_antes = eventos(n).end;
       
       start_depois = eventos(n+1).start;
       end_depois = eventos(n+1).end;
       
       if end_antes >= start_depois
           
           overlap = true;
           
       end
                
    end
    

end

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

end
