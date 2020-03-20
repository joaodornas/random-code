function spikeData = invertSpikeTrain(spike_train,start_time,nBins,bin_size)

tic

disp('BEGIN');

%%%   BUILD BINS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Begin invert spike train...');

disp('...build bins');

    spikes_bins(1:nBins) = struct('bin',[]);

    nSpikes = 0;
    
    spike_train = sort(spike_train);
    
    for j=1:nBins

        spikes_bins(j).bin = sort(spike_train(spike_train>=(((j-1)*bin_size) + start_time) & spike_train<((j*bin_size) + start_time)));

        spikes_bins(j).bin_origin = spikes_bins(j).bin - (j-1)*bin_size - start_time;
        
        nSpikes = nSpikes + length(spikes_bins(j).bin_origin);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   INVERT SPIKE TRAIN   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('...invert spike train');

    if size(spikes_bins,1) == 1
        
        dim = 2;
        
    else
        
        dim = 1;
        
    end
    
    spikes_bins_reversed = flipdim(spikes_bins,dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   RE-BUILD SPIKE TRAIN   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('...re-build spike train');

    spikes = zeros(1,nSpikes);
    
    for j=1:nBins

        spikes_bins_reversed(j).bin = spikes_bins_reversed(j).bin_origin + (j-1)*bin_size + start_time;

        spikes = [spikes spikes_bins_reversed(j).bin]; 

    end

    spikes = spikes(spikes>0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spikeData.spikes_bins = spikes_bins;
spikeData.spikes_bins_reversed = spikes_bins_reversed;
spikeData.spikes = sort(spikes);

disp('END');

toc

end

