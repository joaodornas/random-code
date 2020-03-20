function gallant2002TI(date,site_indeh,channel,registro,start_time,end_time,video_indeh,bin_size,nConditions)


tic

disp('BEGIN');

%%% READ TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Read trials...');

    Spass = load(strcat('_',registro,'-','v',int2str(video_indeh),'.mat'));
    
    for i=1:nConditions
   
        trials_labels(i).label = find(Spass.stimIds == i); 
    
    end
    
    spike_times = Spass.spike_times./ 32000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/',date,'/','sitio',int2str(site_indeh),'/',channel,'/','v',int2str(video_indeh));
    mkdir(filepath,'TI');
    filepath = strcat(filepath,'/','TI','/','v',int2str(video_indeh),'-',date,'-','sitio',int2str(site_indeh),'-',channel);
    filepath = strcat(filepath,'-bin_size-',int2str(bin_size),'-TI');
    

for i=1:nConditions
    
disp('Begin Conditions...');

   trials_spikes = spike_times(trials_labels(i).label,:);
   
   nTrials = size(trials_spikes,1);
   
   nBins = ( end_time - start_time ) / bin_size;
   
   %%%%%  BINARIZA E DIVIDE O SPIKE TRAIN EM WORDS
   
   for t=1:nTrials
       
       spike_train = trials_spikes(t,:);
       spike_train = spike_train(spike_train>0);
       spike_train = spike_train(spike_train>=(start_time/1000) & spike_train<=(end_time/1000));
       
       for b=1:nBins
          
           spikes = size(spike_train(spike_train>=(start_time/1000 + (b-1)*bin_size/1000) & spike_train<(start_time/1000 + (b)*bin_size/1000)),2);
           TI.condicao(i).Trial(t).binWord(b) = spikes;
           
       end
          
   end
   
   for h=1:4
   
       
       disp('Calculate Entropy of Response');
       %%%%% CALCULA AS FREQUENCIAS DAS RESPOSTAS ACROSS TRIALS E A ENTROPIA DA
       %%%%% RESPOSTA

       totalWord = (nTrials/h) * nBins;

       allWords = [];

       for t=1:(nTrials/h)

           for b=1:nBins

                allWords = [allWords, TI.condicao(i).Trial(t).binWord(b)];

           end

       end

       [TI.condicao(i).N(h).WordUnique TI.condicao(i).N(h).nWordUnique] = count_unique(allWords);

       TI.condicao(i).N(h).frequencies = TI.condicao(i).N(h).nWordUnique./totalWord;

       nFrequencies = size(TI.condicao(i).N(h).frequencies,1);

       TI.condicao(i).N(h).HResponse = 0;

       for f=1:nFrequencies

           TI.condicao(i).N(h).HResponse = TI.condicao(i).N(h).HResponse - TI.condicao(i).N(h).frequencies(f)*log2(TI.condicao(i).N(h).frequencies(f));

       end

       sumOfFrequencies = 0;
       for f=1:nFrequencies
          
           sumOfFrequencies = sumOfFrequencies + TI.condicao(i).N(h).frequencies(f)^2;
           
       end
       
       TI.condicao(i).N(h).HResponseMAbound = - log2(sumOfFrequencies);
       
       disp('Calculate Conditional Entropy of Stimulus'); 
       %%%%% CALCULA AS FREQUENCIAS DAS RESPOSTAS POR CONDICAO E A ENTROPIA DA
       %%%%% RESPOSTA
       %%%%% CONDICIONADA AO ESTIMULO


       for b=1:nBins

           totalWord = nTrials/h;

           allWords = [];

           for t=1:(nTrials/h)

                allWords = [allWords, TI.condicao(i).Trial(t).binWord(b)];

           end

           [TI.condicao(i).N(h).Frame(b).WordUnique TI.condicao(i).N(h).Frame(b).nWordUnique] = count_unique(allWords);

           TI.condicao(i).N(h).Frame(b).frequencies = TI.condicao(i).N(h).Frame(b).nWordUnique./totalWord;

           nFrequencies = size(TI.condicao(i).N(h).Frame(b).frequencies,1);

           TI.condicao(i).N(h).Frame(b).HResCond = 0;

           for f=1:nFrequencies

               TI.condicao(i).N(h).Frame(b).HResCond = TI.condicao(i).N(h).Frame(b).HResCond - TI.condicao(i).N(h).Frame(b).frequencies(f)*log2(TI.condicao(i).N(h).Frame(b).frequencies(f));

           end

           sumOfFrequencies = 0;
           for f=1:nFrequencies
          
                sumOfFrequencies = sumOfFrequencies + TI.condicao(i).N(h).Frame(b).frequencies(f)^2;
           
           end
       
           TI.condicao(i).N(h).Frame(b).HResCondMAbound = - log2(sumOfFrequencies);
       
       end

       for b=1:nBins

           HResCond(b) = TI.condicao(i).N(h).Frame(b).HResCond;
           HResCondMAbound(b) = TI.condicao(i).N(h).Frame(b).HResCondMAbound;
           
       end

       disp('Calculate last parameters...');
       
       TI.condicao(i).N(h).meanHResCondMAbound = mean(HResCondMAbound);
       TI.condicao(i).N(h).stdHResCondMAbound = std(HResCondMAbound);
       
       TI.condicao(i).N(h).meanHResCond = mean(HResCond);
       TI.condicao(i).N(h).stdHResCond = std(HResCond);

       TI.condicao(i).N(h).MutualInfo = TI.condicao(i).N(h).HResponse - TI.condicao(i).N(h).meanHResCond;
       TI.condicao(i).N(h).CodingEfficiency = TI.condicao(i).N(h).MutualInfo / TI.condicao(i).N(h).HResponse * 100;
       
       rate = zeros(1,nBins);
       
       for t=1:(nTrials/h)
           
           for b=1:nBins

               rate(b) = rate(b) + TI.condicao(i).Trial(t).binWord(b);

           end

       end
       
       TI.condicao(i).N(h).count = rate;
       TI.condicao(i).N(h).ratePerSec = TI.condicao(i).N(h).count ./ ((nTrials/h)*(bin_size/1000));
       TI.condicao(i).N(h).meanRate = mean(TI.condicao(i).N(h).ratePerSec);
       
       TI.condicao(i).N(h).InfoPerSecond = TI.condicao(i).N(h).MutualInfo / (bin_size/1000) ;
       TI.condicao(i).N(h).InfoPerSpike = TI.condicao(i).N(h).InfoPerSecond / TI.condicao(i).N(h).meanRate ;

   end
   
   
   for h=1:4
      
       HCondEhp(h) = TI.condicao(i).N(h).meanHResCond; 
       HResEhp(h) = TI.condicao(i).N(h).HResponse;
       
   end
   
   f = figure;
   
   plot(1/nTrials/4:1/nTrials/4:1/nTrials,HResEhp,'ro');
   hold on;
   plot(1/nTrials/4:1/nTrials/4:1/nTrials,HCondEhp,'bo');
   hold on;
   
   x = (1/nTrials/4:1/nTrials/4:1/nTrials).';
   s = fitoptions('Method','NonLinearLeastSquares','Lower',[0,0],'Upper',[Inf,Inf],'StartPoint',[0,0,0]);
   g = fittype('Htrue + c1/x + c2/x^2','coeff',{'Htrue','c1','c2'},'options',s);
   length(x)
   length(HResEhp)
   [HRestrue, gofRes] = fit(x,HResEhp.',g);
   [HCondtrue, gofCond] = fit(x,HCondEhp.',g);
   
   plot(HRestrue,'m');
   hold on;
   plot(HCondtrue,'c');
   
   print(f,'-djpeg',strcat(filepath,strcat('-curve-fitting-condition-',int2str(i))));
   
   legend('HResponse','HConditioned');
   
   TI.condicao(i).HRestrue = HRestrue;
   TI.condicao(i).HCondtrue = HCondtrue;
   
   TI.condicao(i).gofRes = gofRes;
   TI.condicao(i).gofCond = gofCond;
   
   TI.condicao(i).HResponse = TI.condicao(i).HRestrue.Htrue;
   TI.condicao(i).HConditionedResponse = TI.condicao(i).HCondtrue.Htrue;
   
   TI.condicao(i).MutualInformation = TI.condicao(i).HResponse - TI.condicao(i).HConditionedResponse;
   
   TI.condicao(i).CodingEfficiency = TI.condicao(i).MutualInformation / TI.condicao(i).HResponse * 100;
    
   TI.condicao(i).InfoPerSecond = TI.condicao(i).MutualInformation / (bin_size/1000);
   
%    rate = 0;
%    for h=1:4
%       
%        rate = rate + TI.condicao(i).N(h).meanRate;
%        
%    end
   
   TI.condicao(i).meanRate = TI.condicao(i).N(1).meanRate;
   
   TI.condicao(i).InfoPerSpike = TI.condicao(i).InfoPerSecond / TI.condicao(i).meanRate;
   
end
       
    disp('Save file.');
    save(filepath,'TI');



end