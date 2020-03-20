function h = simplesEcomplexa(date,site_index,channel,registro,prefFreq)


DTC = load(strcat(registro,'.mat'));

DTC.parameters.stim_duration = 3000;

[MRAllconditions,f1_allConditions,Freq_F1,f0_allConditions,classification,frequency_allConditions,amplitude_allConditions] = simpleComplexClassAllTrialsJP(DTC.psth,DTC.parameters,prefFreq);

spikeCountMatrix = spikecounts(DTC.spike_times,DTC.stimIds,DTC.parameters);

[CTResponse] = cenTendResponse(spikeCountMatrix,'mean');

condition = find(CTResponse==max(CTResponse));

if MRAllconditions(condition) > 1
    
    cellType = 'Simples';
    
else
    
    cellType = 'Complexa';
    
    
end

%filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/');
filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/');

file = fopen(strcat(filepath,cellType,'.txt'),'w');

responseSimplesComplexa = struct('DTC', DTC, 'CTResponse', CTResponse, 'MRAllconditions', MRAllconditions, 'spikeCountMatrix', spikeCountMatrix, 'f1_allConditions', f1_allConditions, 'f0_allConditions', f0_allConditions, 'classification', classification, 'amplitude_allConditions', amplitude_allConditions,'Freq_F1', Freq_F1);

save(strcat(filepath,'responseSimplesComplexa'),'responseSimplesComplexa');


end