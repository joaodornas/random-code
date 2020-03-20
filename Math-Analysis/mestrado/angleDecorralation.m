function decorrelationPairs = angleDecorralation(start_time,end_time)

Spass(1) = load('_nsp015a01_2a-v1.mat');
Spass(2) = load('_nsp015a02_2b-v3.mat');

Spass(3) = load('_nsp016a01_2b-v2.mat');
Spass(4) = load('_nsp016a02_2b-v4.mat');

Spass(5) = load('_nsp017a01_1a-v3.mat');

Spass(6) = load('_nsp017b01_2b-v2.mat');

Spass(7) = load('_nsp018b01_1b-v3.mat');


for k=1:7
    
    Spass(k).spike_times = Spass(k).spike_times./ 32000;
    
    for i=1:2

        trials_label(k,i).label = find(Spass(k).stimIds == i); 
    
        trials_spikes = Spass(k).spike_times(trials_label(k,i).label(:),:);
    
        nTrials(k,i) = size(trials_spikes,1);
    
        
        for T=1:nTrials(k,i)
            
           spikes = trials_spikes(T,:);
                 
           spikes = spikes(spikes>0);
           spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));
           spikes = sort(spikes);
           
           spikes = spikes*1000;
           spikes = round(spikes);
           nSpikes = length(spikes);
                 
           decorrelationPairs.Protocolo(k).RF(i).spike_trains(T).spikes = zeros(1,end_time);
           
           for l=1:nSpikes
           
               decorrelationPairs.Protocolo(k).RF(i).spike_trains(T).spikes(uint8(spikes(l))) = spikes(l);
               
           end
                      
        end
        
    end   

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TEM QUE SER PARES DE NEURONIOS ESTIMULADOS COM O MESMO PROTOCOLO, O MESMO FILME


%%decorrelationaPairs.P13.RF1

nTrialsi = size(decorrelationPairs.Protocolo(1).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(3).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P13.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(1).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(3).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(1).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(3).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P13.RF1.meanAngles = mean(decorrelationPairs.P13.RF1.Angles);


%decorrelationaPairs.P13.RF2

nTrialsi = size(decorrelationPairs.Protocolo(1).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(3).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P13.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(1).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(3).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(1).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(3).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P13.RF2.meanAngles = mean(decorrelationPairs.P13.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P14.RF1

nTrialsi = size(decorrelationPairs.Protocolo(1).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(4).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P14.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(1).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(4).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(1).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(4).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P14.RF1.meanAngles = mean(decorrelationPairs.P14.RF1.Angles);

%decorrelationaPairs.P14.RF2

nTrialsi = size(decorrelationPairs.Protocolo(1).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(4).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P14.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(1).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(4).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(1).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(4).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P14.RF2.meanAngles = mean(decorrelationPairs.P14.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%decorrelationaPairs.P23.RF1

nTrialsi = size(decorrelationPairs.Protocolo(2).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(3).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P23.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(2).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(3).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(2).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(3).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P23.RF1.meanAngles = mean(decorrelationPairs.P23.RF1.Angles);


%decorrelationaPairs.P23.RF2

nTrialsi = size(decorrelationPairs.Protocolo(2).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(3).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P23.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(2).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(3).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(2).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(3).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P23.RF2.meanAngles = mean(decorrelationPairs.P23.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%decorrelationaPairs.P24.RF1

nTrialsi = size(decorrelationPairs.Protocolo(2).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(4).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P24.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(2).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(4).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(2).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(4).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P24.RF1.meanAngles = mean(decorrelationPairs.P24.RF1.Angles);



%decorrelationaPairs.P24.RF2

nTrialsi = size(decorrelationPairs.Protocolo(2).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(4).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P24.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(2).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(4).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(2).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(4).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P24.RF2.meanAngles = mean(decorrelationPairs.P24.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%decorrelationaPairs.P15.RF1


nTrialsi = size(decorrelationPairs.Protocolo(1).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(5).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P15.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(1).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(5).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(1).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(5).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P15.RF1.meanAngles = mean(decorrelationPairs.P15.RF1.Angles);



%decorrelationaPairs.P15.RF2

nTrialsi = size(decorrelationPairs.Protocolo(1).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(5).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P15.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(1).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(5).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(1).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(5).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P15.RF2.meanAngles = mean(decorrelationPairs.P15.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P16.RF1

nTrialsi = size(decorrelationPairs.Protocolo(1).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(6).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P16.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(1).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(6).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(1).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(6).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P16.RF1.meanAngles = mean(decorrelationPairs.P16.RF1.Angles);

%decorrelationaPairs.P16.RF2

nTrialsi = size(decorrelationPairs.Protocolo(1).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(6).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P16.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(1).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(6).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(1).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(6).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P16.RF2.meanAngles = mean(decorrelationPairs.P16.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P17.RF1

nTrialsi = size(decorrelationPairs.Protocolo(1).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(7).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P17.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(1).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(7).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(1).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(7).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P17.RF1.meanAngles = mean(decorrelationPairs.P17.RF1.Angles);

%decorrelationaPairs.P17.RF2

nTrialsi = size(decorrelationPairs.Protocolo(1).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(7).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P17.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(1).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(7).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(1).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(7).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P17.RF2.meanAngles = mean(decorrelationPairs.P17.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P25.RF1 

nTrialsi = size(decorrelationPairs.Protocolo(2).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(5).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P25.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(2).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(5).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(2).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(5).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P25.RF1.meanAngles = mean(decorrelationPairs.P25.RF1.Angles);

%decorrelationaPairs.P25.RF2

nTrialsi = size(decorrelationPairs.Protocolo(2).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(5).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P25.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(2).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(5).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(2).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(5).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P25.RF2.meanAngles = mean(decorrelationPairs.P25.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P26.RF1

nTrialsi = size(decorrelationPairs.Protocolo(2).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(6).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P26.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(2).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(6).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(2).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(6).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P26.RF1.meanAngles = mean(decorrelationPairs.P26.RF1.Angles);

%decorrelationaPairs.P26.RF2

nTrialsi = size(decorrelationPairs.Protocolo(2).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(6).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P26.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(2).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(6).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(2).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(6).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P26.RF2.meanAngles = mean(decorrelationPairs.P26.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P27.RF1

nTrialsi = size(decorrelationPairs.Protocolo(2).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(7).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P27.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(2).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(7).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(2).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(7).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P27.RF1.meanAngles = mean(decorrelationPairs.P27.RF1.Angles);

%decorrelationaPairs.P27.RF2

nTrialsi = size(decorrelationPairs.Protocolo(2).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(7).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P27.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(2).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(7).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(2).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(7).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P27.RF2.meanAngles = mean(decorrelationPairs.P27.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P34.RF1

nTrialsi = size(decorrelationPairs.Protocolo(3).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(4).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P34.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(3).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(4).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(3).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(4).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P34.RF1.meanAngles = mean(decorrelationPairs.P34.RF1.Angles);

%decorrelationaPairs.P34.RF2

nTrialsi = size(decorrelationPairs.Protocolo(3).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(4).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P34.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(3).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(4).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(3).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(4).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P34.RF2.meanAngles = mean(decorrelationPairs.P34.RF2.Angles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%decorrelationaPairs.P35.RF1

nTrialsi = size(decorrelationPairs.Protocolo(3).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(5).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P35.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(3).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(5).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(3).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(5).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P35.RF1.meanAngles = mean(decorrelationPairs.P35.RF1.Angles);

%decorrelationaPairs.P35.RF2

nTrialsi = size(decorrelationPairs.Protocolo(3).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(5).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P35.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(3).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(5).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(3).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(5).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P35.RF2.meanAngles = mean(decorrelationPairs.P35.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P36.RF1

nTrialsi = size(decorrelationPairs.Protocolo(3).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(6).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P36.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(3).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(6).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(3).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(6).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P36.RF1.meanAngles = mean(decorrelationPairs.P36.RF1.Angles);

%decorrelationaPairs.P36.RF2

nTrialsi = size(decorrelationPairs.Protocolo(3).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(6).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P36.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(3).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(6).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(3).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(6).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P36.RF2.meanAngles = mean(decorrelationPairs.P36.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P37.RF1

nTrialsi = size(decorrelationPairs.Protocolo(3).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(7).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P37.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(3).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(7).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(3).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(7).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P37.RF1.meanAngles = mean(decorrelationPairs.P37.RF1.Angles);

%decorrelationaPairs.P37.RF2

nTrialsi = size(decorrelationPairs.Protocolo(3).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(7).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P37.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(3).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(7).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(3).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(7).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P37.RF2.meanAngles = mean(decorrelationPairs.P37.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P45.RF1

nTrialsi = size(decorrelationPairs.Protocolo(4).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(5).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P45.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(4).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(5).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(4).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(5).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P45.RF1.meanAngles = mean(decorrelationPairs.P45.RF1.Angles);

%decorrelationaPairs.P45.RF2

nTrialsi = size(decorrelationPairs.Protocolo(4).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(5).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P45.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(4).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(5).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(4).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(5).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P45.RF2.meanAngles = mean(decorrelationPairs.P45.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P46.RF1

nTrialsi = size(decorrelationPairs.Protocolo(4).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(6).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P46.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(4).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(6).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(4).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(6).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P46.RF1.meanAngles = mean(decorrelationPairs.P46.RF1.Angles);

%decorrelationaPairs.P46.RF2

nTrialsi = size(decorrelationPairs.Protocolo(4).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(6).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P46.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(4).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(6).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(4).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(6).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P46.RF2.meanAngles = mean(decorrelationPairs.P46.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P47.RF1

nTrialsi = size(decorrelationPairs.Protocolo(4).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(7).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P47.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(4).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(7).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(4).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(7).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P47.RF1.meanAngles = mean(decorrelationPairs.P47.RF1.Angles);

%decorrelationaPairs.P47.RF2

nTrialsi = size(decorrelationPairs.Protocolo(4).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(7).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P47.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(4).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(7).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(4).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(7).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P47.RF2.meanAngles = mean(decorrelationPairs.P47.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P56.RF1

nTrialsi = size(decorrelationPairs.Protocolo(5).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(6).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P56.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(5).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(6).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(5).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(6).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P56.RF1.meanAngles = mean(decorrelationPairs.P56.RF1.Angles);

%decorrelationaPairs.P56.RF2

nTrialsi = size(decorrelationPairs.Protocolo(5).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(6).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P56.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(5).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(6).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(5).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(6).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P56.RF2.meanAngles = mean(decorrelationPairs.P56.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P57.RF1

nTrialsi = size(decorrelationPairs.Protocolo(5).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(7).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P57.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(5).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(7).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(5).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(7).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P57.RF1.meanAngles = mean(decorrelationPairs.P57.RF1.Angles);

%decorrelationaPairs.P57.RF2

nTrialsi = size(decorrelationPairs.Protocolo(5).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(7).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P57.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(5).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(7).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(5).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(7).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P57.RF2.meanAngles = mean(decorrelationPairs.P57.RF2.Angles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decorrelationaPairs.P67.RF1

nTrialsi = size(decorrelationPairs.Protocolo(6).RF(1).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(7).RF(1).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P67.RF1.Angles(i) = acos( dot( decorrelationPairs.Protocolo(6).RF(1).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(7).RF(1).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(6).RF(1).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(7).RF(1).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P67.RF1.meanAngles = mean(decorrelationPairs.P67.RF1.Angles);

%decorrelationaPairs.P67.RF2

nTrialsi = size(decorrelationPairs.Protocolo(6).RF(2).spike_trains,1);
nTrialsj = size(decorrelationPairs.Protocolo(7).RF(2).spike_trains,1);

allpairs = allcomb(nTrialsi,nTrialsj);

nAllPairs = size(allpairs,1);

for i=1:nAllPairs
   
    
    decorrelationPairs.P67.RF2.Angles(i) = acos( dot( decorrelationPairs.Protocolo(6).RF(2).spike_trains(allpairs(i,1)).spikes, decorrelationPairs.Protocolo(7).RF(2).spike_trains(allpairs(i,2)).spikes ) /  ( norm(decorrelationPairs.Protocolo(6).RF(2).spike_trains(allpairs(i,1)).spikes) * norm(decorrelationPairs.Protocolo(7).RF(2).spike_trains(allpairs(i,2)).spikes) ) );
    
    
end

decorrelationPairs.P67.RF2.meanAngles = mean(decorrelationPairs.P67.RF2.Angles);


save('decorrelationPairs.mat','decorrelationPairs');

end