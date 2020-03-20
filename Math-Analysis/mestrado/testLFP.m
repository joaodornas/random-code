function testLFP(data,stimIds)


nAllTRials = size(data,1);

for i=1:2
   
    trials_label(i).labels = find(stimIds == i); 
    condition(i).trials = data(trials_label(i).labels,:);
    
end

params.tapers = [3 5];
params.pad = 0;
params.Fs = 1000;
params.fpass = [0.7 157];
params.err = 0;
params.trialave = 1;
movingwin = [1 0.8];

for i=1:2
    
    [S,t,f] = mtspecgramc(condition(i).trials.',movingwin,params);
    figure
    imagesc(t,f,10*log10(abs(S)));
    hold on
    
end

data1 = condition(1).trials;

save('condition1','data1');