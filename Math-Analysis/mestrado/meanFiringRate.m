function meanFiringRate(date,site_index,channel,registro,video_index,start_time,end_time,bin_size,nConditions)

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
       
       meanFR(i) = ( length(spikes_vector) / ( nTrials ) ) / ( (end_time/1000) - (start_time/1000) ) ;
        
       nBins = ( (end_time/1000) - (start_time/1000) ) / (bin_size/1000);
       
        for k=1:nBins
          
                spikes = length(spikes_vector(spikes_vector>=((k-1)*bin_size/1000 + start_time/1000) & spikes_vector<(k*bin_size/1000 + start_time/1000)));
           
                binRate(k) = spikes/((nTrials)*(bin_size/1000));
                
                condition(i).binRate(k) = binRate(k);

        end
       
       minFR(i) = min(binRate);
       maxFR(i) = max(binRate);
       desvio(i) = std(binRate);
      
       
%        capable = @(x) foomean(x);
%        bootStrap(i).ci = bootci(10000,{capable,binRate},'type','bca');
       
       if (i == 2)

            disp('...invert spike train .');
            
            bin_size_resolution = 1;
            
            nBins_new = (end_time - start_time)/bin_size_resolution;
            spikes_inverted = invertSpikeTrain(spikes_vector,start_time,nBins_new,bin_size_resolution);
            spikes_vector_inverted = [];
            spikes_vector_inverted = spikes_inverted.spikes;
            
            spikes_vector_inverted = sort(spikes_vector_inverted);
       
            spikes_vector_inverted = spikes_vector_inverted(spikes_vector_inverted>0);

            spikes_vector_inverted = spikes_vector_inverted(spikes_vector_inverted>(start_time/1000) & spikes_vector_inverted<(end_time/1000));
       
            meanFR(4) = ( length(spikes_vector_inverted) / ( nTrials ) ) / ( (end_time/1000) - (start_time/1000) ) ;
            
            nBins = ( (end_time/1000) - (start_time/1000) ) / (bin_size/1000);
       
            binRate = [];
            
            for k=1:nBins
          
                spikes = length(spikes_vector_inverted(spikes_vector_inverted>=((k-1)*bin_size/1000 + start_time/1000) & spikes_vector_inverted<(k*bin_size/1000 + start_time/1000)));
           
                binRate(k) = spikes/((nTrials)*(bin_size/1000));
                
                condition(4).binRate(k) = binRate(k);

            end
            
            minFR(4) = min(binRate);
            maxFR(4) = max(binRate);
            desvio(4) = std(binRate);
       
%             capable = @(x) foomean(x);
%             bootStrap(i).ci = bootci(10000,{capable,binRate},'type','bca');
       
       end
       
    end

    
hlillieDireto = lillietest(condition(1).binRate)    
hlillieInvertido = lillietest(condition(4).binRate)

[p h] = signrank(condition(1).binRate,condition(4).binRate)

hlillieDireto = lillietest(condition(1).binRate)    
hlillieAE = lillietest(condition(3).binRate)

[p h] = signrank(condition(1).binRate,condition(3).binRate)

hlillieInvertido = lillietest(condition(4).binRate)
hlillieAE = lillietest(condition(3).binRate)    

[p h] = signrank(condition(4).binRate,condition(3).binRate)
   
    
%filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index),'/');
filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index),'/');

mkdir(filepath,'meanFiringRate');
meanFRPath = strcat(filepath,'meanFiringRate/');

for i=1:nConditions+1
    
    file = fopen(strcat(meanFRPath,int2str(meanFR(i)*100),'HZ-',int2str(i),'.txt'),'w');

end

firing = [1 2 3];

sdError = desvio ./ sqrt(nBins) ;

w = figure;
aHand = axes('parent', w);
hold(aHand, 'on');
colors = ['b','r','g'];
for i = 1:(numel(meanFR)-1)
    bar(i, meanFR(i), 'parent', aHand, 'facecolor', colors(i));
end
set(gca, 'XTick', 1:numel(meanFR), 'XTickLabel', {'Filme Direto','Filme Reverso','Atividade Espontanea'},'FontSize',18);
ylabel('{Taxa de Disparo M\''edia}','interpreter','latex','FontSize',24);
errorbar(meanFR,sdError,'k');
print(w,'-depsc',strcat(meanFRPath,'histograma-SP'));

end

% function mm = foomean(x)
%    
%     y = 0;
%         
%     for i=1:length(x)
% 
%         y = y + x(i);
% 
%     end
% 
%     mm = y / length(x);
%     
% end
