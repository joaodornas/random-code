function  singleNeuronAutoCorr(registro,OS)

Spass = load(registro);

nConditions = 3;

start_time = 500;

end_time = 9500;

PSTH_bin_size = 1;

frame_bin_size = 30;

latency_start_time = 500;

latency_end_time = latency_start_time + 300;

for i=1:nConditions
   
    trials_label(i).labels = find(Spass.stimIds == i); 
    
end

spike_times = Spass.spike_times./ 32000;

for i=1:(nConditions + 1)
    
    if i == (nConditions + 1)
        
        j = 2;
        
    else
        
        j = i;
        
    end
        
        
        trials_spikes = spike_times(trials_label(j).labels(:),:); 

        protocol = registro(2:end-4);

        condition = j;

        latency_time = getLatencyForBack(protocol,condition);

        nTrials = size(trials_spikes,1);

        allSpikes = [];

        for k=1:nTrials

            spikes = trials_spikes(k,:);

            a_e = spikes(spikes>=(0/1000) & spikes<(latency_start_time/1000));

            after_ae = spikes(spikes>=(latency_start_time/1000) & spikes<(end_time/1000)) - latency_time;

            spikes = [];

            spikes = [a_e, after_ae];

            spikes = spikes(spikes>0);

            spikes = sort(spikes);

            spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));

            if (j == (nConditions + 1))

                bin_size = 1;

                nBins = (end_time - start_time)/bin_size;

                %disp(strcat('...invert trial .',int2str(k)));

                spikes_inverted = invertSpikeTrain(spikes,start_time,nBins,bin_size);
                spikes = spikes_inverted.spikes;

            end

            Spike_Trains_Conditions(i).allTrials(k).trial = spikes;

            allSpikes = [allSpikes, spikes];

        end

        allSpikes = reshape(allSpikes.',[],1);

        allSpikes = sort(allSpikes);

        allSpikes = allSpikes(allSpikes>0);

        nBins_PSTH = (end_time - start_time) / PSTH_bin_size;

        nBins_frames = (end_time - start_time) / frame_bin_size;

        for k=1:nBins_PSTH

            PSTH_Conditions_1ms(i).rate(k) = numel(allSpikes(allSpikes>=((start_time/1000) + (k-1)*(PSTH_bin_size/1000)) & allSpikes<=((start_time/1000) + (k)*(PSTH_bin_size/1000)) ) ) / ( nTrials * PSTH_bin_size/1000 );

            PSTH_Conditions_1ms(i).spikes(k).times = allSpikes(allSpikes>=((start_time/1000) + (k-1)*(PSTH_bin_size/1000)) & allSpikes<=((start_time/1000) + (k)*(PSTH_bin_size/1000)) ) ;

            PSTH_Conditions_1ms(i).sd(k) = std(PSTH_Conditions_1ms(i).spikes(k).times);

        end

        for k=1:nBins_frames

            PSTH_Conditions_frames(i).rate(k) = numel(allSpikes(allSpikes>=((start_time/1000) + (k-1)*(frame_bin_size/1000)) & allSpikes<=((start_time/1000) + (k)*(frame_bin_size/1000)) ) ) / ( nTrials * frame_bin_size/1000 );

            PSTH_Conditions_frames(i).spikes(k).times = allSpikes(allSpikes>=((start_time/1000) + (k-1)*(frame_bin_size/1000)) & allSpikes<=((start_time/1000) + (k)*(frame_bin_size/1000)) );

            PSTH_Conditions_frames(i).sd(k) = std(PSTH_Conditions_frames(i).spikes(k).times);

        end
    
end

Correlations_PSTH_1ms_Conditions = getAutoCorrPSTH(PSTH_Conditions_1ms,'PSTH_1ms');
Correlations_PSTH_frames_Conditions = getAutoCorrPSTH(PSTH_Conditions_frames,'PSTH_frames');

Fits_Gauss_1ms_Conditions = getGaussianFitsAutoCorr(Correlations_PSTH_1ms_Conditions,'PSTH_1ms');
Fits_Gauss_frames_Conditions = getGaussianFitsAutoCorr(Correlations_PSTH_frames_Conditions,'PSTH_frames');

Correlations_Spike_Conditions = getAutoCorrSpike(Spike_Trains_Conditions,start_time,end_time);
Fits_Gauss_Spike_Conditions = getGaussianFitsAutoCorr(Correlations_Spike_Conditions,'Spike');

name = registro(2:end-4);

if strcmp(OS,'Win')
    
    folderpath = strcat('Z:\DATA\Forward-Backward\Alonso\singleNeuron\');
    bar = '\';
    
elseif strcmp(OS,'Mac')
    
    folderpath = strcat('/Volumes/Data/DATA/Forward-Backward/Alonso/singleNeuron/');
    bar = '/';
    
end

mkdir(folderpath,name);

results = struct('name',name,'Spike_Trains_Conditions',Spike_Trains_Conditions,'PSTH_Conditions_1ms',PSTH_Conditions_1ms,'PSTH_Conditions_frames',PSTH_Conditions_frames,'Correlations_PSTH_1ms_Conditions',Correlations_PSTH_1ms_Conditions,'Correlations_PSTH_frames_Conditions',Correlations_PSTH_frames_Conditions,'Correlations_Spike_Conditions',Correlations_Spike_Conditions,'Fits_Gauss_1ms_Conditions',Fits_Gauss_1ms_Conditions,'Fits_Gauss_frames_Conditions',Fits_Gauss_frames_Conditions,'Fits_Gauss_Spike_Conditions',Fits_Gauss_Spike_Conditions);

save(strcat(folderpath,name,bar,name,'-','singleNeuronAutoCorr'),'results');

for c=1:(nConditions + 1)
    
    plotFigureFitCorrelation(Fits_Gauss_1ms_Conditions,name,'PSTH_1ms',c,folderpath,bar);

    plotFigureFitCorrelation(Fits_Gauss_frames_Conditions,name,'PSTH_frames',c,folderpath,bar);

    plotFigureFitCorrelation(Fits_Gauss_Spike_Conditions,name,'Spike',c,folderpath,bar);
    
end


function plotFigureFitCorrelation(Fits_Gauss_Conditions,name,kind,c,folderpath,bar)


    if strcmp(kind,'PSTH_1ms')
        
        auto_correlogram = Fits_Gauss_Conditions(c).correlogram_PSTH_1ms.correlogram_norm;
        auto_xcorr = Fits_Gauss_Conditions(c).xcorr_PSTH_1ms.xcorr_norm;
        
        trace = (-length(auto_correlogram)/2 + 1):(length(auto_correlogram)/2);
        tracex = (-length(auto_xcorr)/2 + 1):(length(auto_xcorr)/2);

        A = Fits_Gauss_Conditions(c).correlogram_PSTH_1ms.A;
        Ax = Fits_Gauss_Conditions(c).xcorr_PSTH_1ms.Ax;

        sigma = Fits_Gauss_Conditions(c).correlogram_PSTH_1ms.sigma;
        sigmax = Fits_Gauss_Conditions(c).xcorr_PSTH_1ms.sigmax;

        mu = Fits_Gauss_Conditions(c).correlogram_PSTH_1ms.mu;
        mux = Fits_Gauss_Conditions(c).xcorr_PSTH_1ms.mux;
        
    elseif strcmp(kind,'PSTH_frames')
        
        auto_correlogram = Fits_Gauss_Conditions(c).correlogram_PSTH_frames.correlogram_norm;
        auto_xcorr = Fits_Gauss_Conditions(c).xcorr_PSTH_frames.xcorr_norm;
        
        trace = (-length(auto_correlogram)/2 + 1):(length(auto_correlogram)/2);
        tracex = (-length(auto_xcorr)/2 + 1):(length(auto_xcorr)/2);

        A = Fits_Gauss_Conditions(c).correlogram_PSTH_frames.A;
        Ax = Fits_Gauss_Conditions(c).xcorr_PSTH_frames.Ax;

        sigma = Fits_Gauss_Conditions(c).correlogram_PSTH_frames.sigma;
        sigmax = Fits_Gauss_Conditions(c).xcorr_PSTH_frames.sigmax;

        mu = Fits_Gauss_Conditions(c).correlogram_PSTH_frames.mu;
        mux = Fits_Gauss_Conditions(c).xcorr_PSTH_frames.mux;
        
    elseif strcmp(kind,'Spike')
        
        auto_correlogram = Fits_Gauss_Conditions(c).correlogram_Spike.correlogram_norm;
        auto_xcorr = Fits_Gauss_Conditions(c).xcorr_Spike.xcorr_norm;
        
        trace = (-length(auto_correlogram)/2 + 1):(length(auto_correlogram)/2);
        tracex = (-length(auto_xcorr)/2 + 1):(length(auto_xcorr)/2);

        A = Fits_Gauss_Conditions(c).correlogram_Spike.A;
        Ax = Fits_Gauss_Conditions(c).xcorr_Spike.Ax;

        sigma = Fits_Gauss_Conditions(c).correlogram_Spike.sigma;
        sigmax = Fits_Gauss_Conditions(c).xcorr_Spike.sigmax;

        mu = Fits_Gauss_Conditions(c).correlogram_Spike.mu;
        mux = Fits_Gauss_Conditions(c).xcorr_Spike.mux;
        
    end
        

        f = figure;

        x = trace;
        xx = tracex;

        y = A * exp((-(x - mu).^2)/(2*sigma^2));
        yy = Ax * exp((-(xx - mux).^2)/(2*sigmax^2));

        y = y ./ max(y);
        yy = yy ./ max(yy);

        plot(trace,auto_correlogram,'r')
        hold on;
        plot(x,y,'b');
        text(0.7*max(trace),0.9*max(auto_correlogram),strcat('SIGMA:',num2str(sigma)),'FontSize',20);
        text(0.7*max(trace),0.7*max(auto_correlogram),strcat('mu:',num2str(mu)),'FontSize',20);
        print(f,'-depsc',strcat(folderpath,name,bar,name,'-','singleNeuronAutoCorr-c',int2str(c),'-sigma-',num2str(sigma),'-mu-',num2str(mu),'-',kind,'.eps'));
        close all;

        g = figure;

        plot(tracex,auto_xcorr,'r')
        hold on;
        plot(xx,yy,'b');
        text(0.7*max(tracex),0.9*max(auto_xcorr),strcat('SIGMA:',num2str(sigmax)),'FontSize',20);
        text(0.7*max(tracex),0.7*max(auto_xcorr),strcat('mu:',num2str(mux)),'FontSize',20);
        print(g,'-depsc',strcat(folderpath,name,bar,name,'-','singleNeuronAutoCorr-c',int2str(c),'-sigmax-',num2str(sigmax),'-mux-',num2str(mux),'-',kind,'.eps'));
        close all;
        
end
    
end





