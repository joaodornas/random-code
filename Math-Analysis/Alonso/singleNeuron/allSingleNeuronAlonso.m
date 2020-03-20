function allSingleNeuronAlonso(OS)

registro = importdata('memoryBackwardProtocols.txt');

for r=1:length(registro)

    singleNeuronAutoCorr(registro{r},OS);

end

end

