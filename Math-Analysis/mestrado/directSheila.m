function directSheila(date,site_index,channel,registro,video_index,start_time,end_time,resolution,bin_size)

tic

disp('BEGIN');

%%% READ TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Read trials...');

    Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));
    
    for i=1:3
   
        trials_label(i).label = find(Spass.stimIds == i); 
    
    end
    
    Spass.spike_times = Spass.spike_times./ 32000;
    
    for i=1:3
   
        memoryFrameData(i).trials_spikes = Spass.spike_times(trials_label(i).label,:);

        nTrials(i) = size(memoryFrameData(i).trials_spikes,1);

        for k=1:nTrials(i)

            spikes = memoryFrameData(i).trials_spikes(k,:);
            spikes = spikes(spikes>0);
            spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));

            memoryFrameData(i).trials(k).spikes = sort(spikes);
            
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CALCULATE THE NUMBER OF WORDS E LETTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculate the number of words and letters...');

%... letters
nLetters = bin_size/resolution;
lLetters = resolution/1000;

for i=1:3
    
    memoryFrameData(i).nLetters = nLetters;
    memoryFrameData(i).lLetters = lLetters;

end

%... words
nWords = (end_time - start_time)/bin_size;
lWords = bin_size/1000;

for i=1:3
    
    memoryFrameData(i).nWords = nWords;
    memoryFrameData(i).lWords = lWords;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DIVIDE IN BINS/WORDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Devide in Bins/Words...');

for i=1:3

    for k=1:nTrials(i)
        
        spikes = memoryFrameData(i).trials(k).spikes;
        
        for w=1:nWords
            
            for l=1:nLetters
                
                nSpikes(l) = size(spikes(spikes>=(start_time/1000 + (w-1)*lWords + (l-1)*lLetters) & spikes<(start_time/1000 + (w-1)*lWords + l*lLetters)),2);
        
            end
            
            word = 0;
            
            for l=nLetters:-1:1
               
                word = word + nSpikes(l)*10^(nLetters-l);
                
            end
            
             memoryFrameData(i).trials_words(k).word(w) = word;
            
            clear nSpikes word;
            
        end
        
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l=0:(nWords - 1)
    
    newNWords = nWords - l;
    
%%% CALCULA A FREQUENCIA DE CADA WORD DENTRO DE TODOS OS TRIALS %%%%%%%%%%%
disp('Calculate the frequency of each word across all trials...');

    for i=1:3

        totalWords = nTrials(i)*newNWords;

        allWords = [];

        for k=1:nTrials(i)

            words =  memoryFrameData(i).trials_words(k).word(:);
            
            if l~=0
                
                words(1:l) = [];
            
            end
            
            allWords = [allWords, words];

        end

        [memoryFrameData(i).l(l+1).unique  memoryFrameData(i).l(l+1).nUnique] = count_unique(allWords);

        clear allWords;

        memoryFrameData(i).l(l+1).frequencies = memoryFrameData(i).l(l+1).nUnique./totalWords;

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% CALCULA A FREQUENCIA DE CADA WORD DENTRO DE TODOS OS TRIALS %%%%%%%%%%%
%%% CONDICIONADO A CADA FRAME/ESTIMULO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculate each word frequency conditioned to a frame...');

    for i=1:3

        for w=(1+l):nWords

            condWords = [];

            for k=1:nTrials(i)

                condWords = [condWords, memoryFrameData(i).trials_words(k).word(w)];

            end

            [memoryFrameData(i).l(l+1).wFrames(w).unique memoryFrameData(i).l(l+1).wFrames(w).nUnique] = count_unique(condWords);

            memoryFrameData(i).l(l+1).wFrames(w).frequencies = memoryFrameData(i).l(l+1).wFrames(w).nUnique./nTrials(i);

            clear condWords;

        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CALCULA A ENTROPIA DA RESPOSTA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculate the response entropy...');

    for i=1:3

        memoryFrameData(i).l(l+1).HResponse = 0;

        nFrequencies = size(memoryFrameData(i).l(l+1).frequencies,1);

        for f=1:nFrequencies

            memoryFrameData(i).l(l+1).HResponse = memoryFrameData(i).l(l+1).HResponse - memoryFrameData(i).l(l+1).frequencies(f)*log2(memoryFrameData(i).l(l+1).frequencies(f));

        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CALCULA A ENTROPIA DA RESPOSTA CONDICIONADA AO ESTIMULO %%%%%%%%%%%%%%%
disp('Calculate the response entropy conditioned to a frame...');

    for i=1:3

        memoryFrameData(i).l(l+1).HCondResponse = 0;

        for w=(l+1):nWords

            memoryFrameData(i).l(l+1).wFrames(w).HCondResponse = 0;

            nFrequencies = size(memoryFrameData(i).l(l+1).wFrames(w).frequencies,2);

            for f=1:nFrequencies

                memoryFrameData(i).l(l+1).wFrames(w).HCondResponse = memoryFrameData(i).l(l+1).wFrames(w).HCondResponse - memoryFrameData(i).l(l+1).wFrames(w).frequencies(f)*log2(memoryFrameData(i).l(l+1).wFrames(w).frequencies(f));

            end

        end

        for w=(l+1):nWords

            HCondResponse(w) = memoryFrameData(i).l(l+1).wFrames(w).HCondResponse;

        end

        memoryFrameData(i).l(l+1).HCondResponse = mean(HCondResponse,2);
        memoryFrameData(i).l(l+1).HCondResponseStd = std(HCondResponse,0,2);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CALCULA INFORMACAO MUTUA E EFICIENCIA DE CODIFICACAO %%%%%%%%%%%%%%%%%%
disp('Calculate mutual information and coding efficiency...');

    for i=1:3

        memoryFrameData(i).l(l+1).MutualInfo = memoryFrameData(i).l(l+1).HResponse - memoryFrameData(i).l(l+1).HCondResponse;

        memoryFrameData(i).l(l+1).CodingEfficiency = memoryFrameData(i).l(l+1).MutualInfo / memoryFrameData(i).l(l+1).HResponse;
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   CALCULATE THE MAXIMUM CODING EFFICIENCY   %%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculate the maximum coding efficiency...');

for i=1:3

    allCodingEfficiency = [];

    for l=0:(nWords-1)
    
        allCodingEfficiency = [allCodingEfficiency, memoryFrameData(i).l(l+1).CodingEfficiency];    
    
    end
    
    memoryFrameData(i).maxCodingEfficiency = max(allCodingEfficiency);
    memoryFrameData(i).maxL = find(allCodingEfficiency==memoryFrameData(i).maxCodingEfficiency);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% SALVA NO DISCO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save data...');

filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));

filepath = strcat(filepath,'/','memoria-frame','/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel);

filepath = strcat(filepath,'-res-',int2str(resolution),'-bin-',int2str(bin_size),'-memoria-frame');
    
save(filepath,'memoryFrameData');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% GERA GRAFICO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plot graphs...');

for i=1:3

    CodingEfficiency = [];
    
    for l=0:(nWords-1)
       
        CodingEfficiency = [CodingEfficiency, memoryFrameData(i).l(l+1).CodingEfficiency];    
        
    end
    
    f(i) = figure;
    
    bar(CodingEfficiency);
    
    print(f(i),'-dbmp',strcat(filepath,'-bar-condition-',int2str(i),'.bmp'));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

disp('END');

toc

end
