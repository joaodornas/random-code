function vonMiseFourier

registro = '_nsp008a03_2a-v4.mat';



[X, Y, R] = receptiveField(registro);

v = 4;

name = registro;

nFrames = 300;

[ForwardMovieCRF, BackwardMovieCRF] = ForBackCRFmovie(name,v,300,X,Y,68,0);





Spass = load(char(registro));

forwardlabels = find(Spass.stimIds == 1); 

backwardlabels = find(Spass.stimIds == 2); 

spike_times = Spass.spike_times ./ 32000;

forwardtrials = spike_times(forwardlabels(:),:);

backwardtrials = spike_times(backwardlabels(:),:);

start_time = 500;

end_time = 9500;

bin_size = 30;

PSTH_validation_FOR = getPSTH_validation(forwardtrials,start_time,end_time,bin_size);

PSTH_validation_BACK = getPSTH_validation(backwardtrials,start_time,end_time,bin_size);





DTC_results = load('/Volumes/Data/DATA/Forward-Backward/DTC/preferida.mat');

preferida = DTC_results.results(1,17).preferida;
anti_preferida = DTC_results.results(1,17).anti_preferida;
DI = DTC_results.results(1,17).DI;
A = DTC_results.results(1,17).A;
m = DTC_results.results(1,17).m;





for g=1:nFrames
    
   For_fourier(g).cdata = fft2(ForwardMovieCRF(g).cdata);
   Back_fourier(g).cdata = fft2(BackwardMovieCRF(g).cdata);
   
   Ang_For(g).cdata = angle(For_fourier(g).cdata);
   Ang_Back(g).cdata = angle(Back_fourier(g).cdata);
        
end

Dir_For(1).cdata = Ang_For(1).cdata;
Dir_Back(1).cdata = Ang_Back(1).cdata;

for g=2:nFrames
   
    Dir_For(g).cdata = Ang_For(g).cdata - Ang_For(g-1).cdata;
    
    Dir_Back(g).cdata = Ang_For(g).cdata - Ang_For(g-1).cdata;
    
end


PSTH_prediction_FOR = getPSTH_prediction(Dir_For,preferida,anti_preferida,DI,A,nFrames);

PSTH_prediction_BACK = getPSTH_prediction(Dir_Back,preferida,anti_preferida,DI,A,nFrames);


[r_for , p_for] = corr(PSTH_prediction_FOR,PSTH_validation_FOR,'type','pearson');

[r_back , p_back] = corr(PSTH_prediction_BACK,PSTH_validation_BACK,'type','pearson');

%PSTH_validation_FOR = PSTH_validation_FOR ./ max(PSTH_validation_FOR);

%PSTH_prediction_FOR = PSTH_prediction_FOR ./ max(PSTH_prediction_FOR);


f = figure;
plot(1:nFrames,PSTH_validation_FOR,'b');
hold on;
plot(1:nFrames,PSTH_prediction_FOR,'r');
text(nFrames,max(PSTH_validation_FOR),num2str(r_for));
print(f,'-depsc',strcat('/Volumes/Data/_Research/_doutorado/models/PSTH-validation-prediction-Forward'));

%PSTH_validation_BACK = PSTH_validation_BACK ./ max(PSTH_validation_BACK);

%PSTH_prediction_BACK = PSTH_prediction_BACK ./ max(PSTH_prediction_BACK);

g = figure;
plot(1:nFrames,PSTH_validation_BACK,'b');
hold on;
plot(1:nFrames,PSTH_prediction_BACK,'r');
text(nFrames,max(PSTH_validation_BACK),num2str(r_back));
print(f,'-depsc',strcat('/Volumes/Data/_Research/_doutorado/models/PSTH-validation-prediction-Backward'));


function PSTH = getPSTH_prediction(movie,preferida,anti_preferida,DI,A,nFrames)

    MM = size(movie(1).cdata,1);
    NN = size(movie(1).cdata,2);
    
    k = 2;
    k_2 = 2;

    if (DI < 0.5) || (DI > 0.9)

        A_2 = 0;

    else

        A_2 = 0.7 * A;

    end

    
    
    for n=1:nFrames
        
        M = vonMiseMatrix(movie(n).cdata,A,preferida,anti_preferida,DI);
        
        PSTH(n) = M;
        
    end
    

end

function PSTH = getPSTH_validation(trials,start_time,end_time,bin_size)
    
    nTrials = size(trials,1);
        
    allTrials = reshape(trials.',[],1);
    
    allTrials = allTrials(allTrials>0);
    
    allTrials = sort(allTrials);
    
    nBins = (end_time - start_time)/bin_size;
    
    for k=1:nBins

        spikes = numel(allTrials(allTrials>=((start_time/1000) + (((k-1)*bin_size/1000))) & allTrials<((start_time/1000) + k*(bin_size/1000))));
        
        PSTH(k) = spikes / ( nTrials * (bin_size/1000) ) ;
        
    end

end

end

