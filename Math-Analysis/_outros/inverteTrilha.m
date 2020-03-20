function spikes = inverteTrilha(spike_train,start_time,nBins)

%disp('BEGIN - inverteTrilha');

%%%  DEFINE TAMANHO DO BIN   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_size = 1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  CONSTROI UMA TRILHA DE BINS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
spikes_bins(1:nBins) = struct('bin',[]);

nSpikes = 0;

spike_train = spike_train(spike_train>0);

spike_train = sort(spike_train);

for j=1:nBins

    spikes_bins(j).bin = sort(spike_train(spike_train>=(((j-1)*bin_size/1000) + start_time/1000) & spike_train<((j*bin_size/1000) + start_time/1000)));

    spikes_bins(j).bin_origin = spikes_bins(j).bin - (j-1)*bin_size/1000 - start_time/1000;

    nSpikes = nSpikes + length(spikes_bins(j).bin_origin);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   INVERTE TRILHA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(spikes_bins,1) == 1

    dim = 2;

else

    dim = 1;

end

spikes_bins_reversed = flipdim(spikes_bins,dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   RECONSTROI A TRILHA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spikes = zeros(1,nSpikes);

for j=1:nBins

    spikes_bins_reversed(j).bin = spikes_bins_reversed(j).bin_origin + (j-1)*bin_size/1000 + start_time/1000;

    spikes = [spikes spikes_bins_reversed(j).bin]; 

end

spikes = spikes(spikes>0);

spikes = sort(spikes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp('END - inverteTrilha');

end

