
%%% get only PRIMATES PMIDs

load('__PMIDs-Retina_v4b.mat');
fascicularis = animal.macaca_fascicularis.PMID;
mulatta = animal.macaca_mulatta.PMID;
nemestrina = animal.macaca_nemestrina.PMID;
rhesus = animal.rhesus.PMID;
primate = animal.primate.PMID;
monkey = animal.monkey.PMID;
macaque = animal.macaque.PMID;
human = animal.human.PMID;
all_primates = [fascicularis, mulatta, nemestrina, rhesus, primate, monkey, macaque, human];
all_primates = unique(all_primates);

save('all_primates.mat','all_primates');

primates_labels = {'macaca_fascicularis', 'macaca_mulatta', 'macaca_nemestrina', 'rhesus', 'primate', 'monkey', 'macaque', 'human'};
all_PMID = [];
for iP=1:length(primates_labels)

    eval(strcat('these_fields = fields(animal_neuron.',primates_labels{iP},');'));

    for iF=1:length(these_fields)
    
        eval(strcat('PMID = animal_neuron.',primates_labels{iP},'.',these_fields{iF},'.PMID;'));

        all_PMID = [all_PMID, PMID];
    
    end

end

all_PMID = unique(all_PMID);

save('all_primates_neuron.mat','all_PMID');
