function corrigirRegistro(registro,video_index,spikes)

Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));

Spass.stimIds(spikes) = [];
Spass.spike_times(spikes,:) = [];

save(strcat('_',registro,'-','v',int2str(video_index),'-corrigido','.mat'),'-struct','Spass');

end

