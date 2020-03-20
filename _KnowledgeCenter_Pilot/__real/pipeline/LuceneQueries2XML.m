
%%% LUCENE QUERIES 2 XML/OPML MINDLY MAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DATABASE
goal = 'Retina';
version = strcat(int2str(5));

load(strcat('__descriptives-',goal,'_v',version,'.mat'));

ANIMAL = animals;
clear animals
PHENOMENON = phenomenon;
clear phenomenon
NEURON = neuron;
clear neuron

nAnimals = length(ANIMAL);
nPhenomenon = length(PHENOMENON);
nNeuron = length(NEURON);

load('query-lucene-Retina-v5-PMID.mat');

nQueries = length(all_queries);

for iQ=1:nQueries
   
    labels{iQ} = all_queries(iQ).label;
    
end

allPMID = [];

% only Animals
for iAnimals=1:nAnimals
    
    disp(strcat('query-lucene-animals-',ANIMAL{iAnimals},'.txt'));

    %%% data = loadjson(strcat('query-lucene-animals-',ANIMAL{iAnimals},'.txt'));
    
    idx = strmatch(strcat('query-lucene-animals-',ANIMAL{iAnimals}),labels,'exact');
    %%% Index = find(not(cellfun('isempty', idx)));
    
    nCount = all_queries(idx).result.response.numFound;
    
    if nCount > 0
        
        eval(strcat('animal.',strrep(ANIMAL{iAnimals},' ','_'),'.count = ',int2str(nCount),';'));
        
        eval(strcat('all.animal.',strrep(ANIMAL{iAnimals},' ','_'),'= ',int2str(nCount),';'));
        
        tmp_PMID = all_queries(idx).result.facet_counts.facet_fields.PMID;
        
        nPMIDs = length(tmp_PMID);
        count = 0;
        for iPMID=1:2:nPMIDs
            count = count + 1;
            PMIDs{count} = tmp_PMID{iPMID};
        end
        
        allPMID = [allPMID,PMIDs];
        
        eval(strcat('animal.',strrep(ANIMAL{iAnimals},' ','_'),'.PMID = ','PMIDs',';'));
        
        clear PMIDs tmp_PMID
        
    end
    
end

% only Phenomenon
for iPhenomenon=1:nPhenomenon
    
    disp(strcat('query-lucene-phenomenon-',PHENOMENON{iPhenomenon},'.txt'));

    %%% data = loadjson(strcat('query-lucene-phenomenon-',PHENOMENON{iPhenomenon},'.txt'));
    
    idx = strmatch(strcat('query-lucene-phenomenon-',PHENOMENON{iPhenomenon}),labels,'exact');
    %%% Index = find(not(cellfun('isempty', idx)));
    
    nCount = all_queries(idx).result.response.numFound;
    
    if nCount > 0

        eval(strcat('phenomenon.',PHENOMENON{iPhenomenon},'.count = ',int2str(nCount),';'));
        
        eval(strcat('all.phenomenon.',PHENOMENON{iPhenomenon},'=',int2str(nCount),';'));
        
        tmp_PMID = all_queries(idx).result.facet_counts.facet_fields.PMID;
        
        nPMIDs = length(tmp_PMID);
        count = 0;
        for iPMID=1:2:nPMIDs
            count = count + 1;
            PMIDs{count} = tmp_PMID{iPMID};
        end
        
        allPMID = [allPMID,PMIDs];
        
        eval(strcat('phenomenon.',PHENOMENON{iPhenomenon},'.PMID = ','PMIDs',';'));
        
        clear PMIDs tmp_PMID
        
    end
    
end

% only Neuron
for iNeuron=1:nNeuron
    
    disp(strcat('query-lucene-neuron-',NEURON{iNeuron},'.txt'));

    %%% data = loadjson(strcat('query-lucene-neuron-',NEURON{iNeuron},'.txt'));
    
    idx = strmatch(strcat('query-lucene-neuron-',NEURON{iNeuron}),labels,'exact');
    %%% Index = find(not(cellfun('isempty', idx)));
    
    nCount = all_queries(idx).result.response.numFound;
    
    if nCount > 0
        
        eval(strcat('neuron.',NEURON{iNeuron},'.count = ',int2str(nCount),';'));
        
        eval(strcat('all.neuron.',NEURON{iNeuron},'=',int2str(nCount),';'));
        
        tmp_PMID = all_queries(idx).result.facet_counts.facet_fields.PMID;
        
        nPMIDs = length(tmp_PMID);
        count = 0;
        for iPMID=1:2:nPMIDs
            count = count + 1;
            PMIDs{count} = tmp_PMID{iPMID};
        end
        
        allPMID = [allPMID,PMIDs];
        
        eval(strcat('neuron.',NEURON{iNeuron},'.PMID = ','PMIDs',';'));
        
        clear PMIDs tmp_PMID
        
    end

end

% Animal & Neuron
for iAnimals=1:nAnimals
    
    for iNeuron=1:nNeuron
               
        disp(strcat('query-lucene-animals-neuron-',strrep(ANIMAL{iAnimals},' ','_'),'-',NEURON{iNeuron},'.txt'));

        %%% data = loadjson(strcat('query-lucene-animals-neuron-',ANIMAL{iAnimals},'-',NEURON{iNeuron},'.txt'));
        
        idx = strmatch(strcat('query-lucene-animals-neuron-',ANIMAL{iAnimals},'-',NEURON{iNeuron}),labels,'exact');
        %%% Index = find(not(cellfun('isempty', idx)));
        
        nCount = all_queries(idx).result.response.numFound;
        
        if nCount > 0
        
            eval(strcat('animal_neuron.',strrep(ANIMAL{iAnimals},' ','_'),'.',NEURON{iNeuron},'.count = ',int2str(nCount),';'));
            
            eval(strcat('all.animal_neuron.',strrep(ANIMAL{iAnimals},' ','_'),'.',NEURON{iNeuron},'=',int2str(nCount),';'));
            
            tmp_PMID = all_queries(idx).result.facet_counts.facet_fields.PMID;
        
            nPMIDs = length(tmp_PMID);
            count = 0;
            for iPMID=1:2:nPMIDs
                count = count + 1;
                PMIDs{count} = tmp_PMID{iPMID};
            end
        
            allPMID = [allPMID,PMIDs];
            
            eval(strcat('animal_neuron.',strrep(ANIMAL{iAnimals},' ','_'),'.',NEURON{iNeuron},'.PMID = ','PMIDs',';'));
            
            clear PMIDs tmp_PMID
            
        end
        
    end

end

% Animal & Neuron & Phenomenon
for iAnimals=1:nAnimals
    
    for iNeuron=1:nNeuron
        
        for iPhenomenon=1:nPhenomenon
            
            disp(strcat('query-lucene-animals-neuron-phenomenon-',strrep(ANIMAL{iAnimals},' ','_'),'-',NEURON{iNeuron},'-',PHENOMENON{iPhenomenon},'.txt'));

            %%% data = loadjson(strcat('query-lucene-animals-neuron-phenomenon-',ANIMAL{iAnimals},'-',NEURON{iNeuron},'-',PHENOMENON{iPhenomenon},'.txt'));
            
            idx = strmatch(strcat('query-lucene-animals-neuron-phenomenon-',ANIMAL{iAnimals},'-',NEURON{iNeuron},'-',PHENOMENON{iPhenomenon}),labels,'exact');
            %%% Index = find(not(cellfun('isempty', idx)));
        
            nCount = all_queries(idx).result.response.numFound;
            
            if nCount > 0
        
                eval(strcat('animal_neuron_phenomenon.',strrep(ANIMAL{iAnimals},' ','_'),'.',NEURON{iNeuron},'.',PHENOMENON{iPhenomenon},'.count = ',int2str(nCount),';'));
                
                eval(strcat('all.animal_neuron_phenomenon.',strrep(ANIMAL{iAnimals},' ','_'),'.',NEURON{iNeuron},'.',PHENOMENON{iPhenomenon},'=',int2str(nCount),';'));
                
                tmp_PMID = all_queries(idx).result.facet_counts.facet_fields.PMID;
        
                nPMIDs = length(tmp_PMID);
                count = 0;
                for iPMID=1:2:nPMIDs
                    count = count + 1;
                    PMIDs{count} = tmp_PMID{iPMID};
                end
                
                allPMID = [allPMID,PMIDs];
                
                eval(strcat('animal_neuron_phenomenon.',strrep(ANIMAL{iAnimals},' ','_'),'.',NEURON{iNeuron},'.',PHENOMENON{iPhenomenon},'.PMID = ','PMIDs',';'));
                
                clear PMIDs tmp_PMID
                
            end
            
        end
        
    end

end

% Neuron & Animal
for iNeuron=1:nNeuron
    
    for iAnimals=1:nAnimals
               
        disp(strcat('query-lucene-neuron-animals-',strrep(ANIMAL{iAnimals},' ','_'),'-',NEURON{iNeuron},'.txt'));

        %%% data = loadjson(strcat('query-lucene-animals-neuron-',ANIMAL{iAnimals},'-',NEURON{iNeuron},'.txt'));
        
        idx = strmatch(strcat('query-lucene-animals-neuron-',ANIMAL{iAnimals},'-',NEURON{iNeuron}),labels,'exact');
        %%% Index = find(not(cellfun('isempty', idx)));
        
        nCount = all_queries(idx).result.response.numFound;
        
        if nCount > 0
        
            eval(strcat('neuron_animal.',NEURON{iNeuron},'.',strrep(ANIMAL{iAnimals},' ','_'),'.count = ',int2str(nCount),';'));
            
            eval(strcat('all.neuron_animal.',NEURON{iNeuron},'.',strrep(ANIMAL{iAnimals},' ','_'),'=',int2str(nCount),';'));
            
            tmp_PMID = all_queries(idx).result.facet_counts.facet_fields.PMID;
        
            nPMIDs = length(tmp_PMID);
            count = 0;
            for iPMID=1:2:nPMIDs
                count = count + 1;
                PMIDs{count} = tmp_PMID{iPMID};
            end
        
            allPMID = [allPMID,PMIDs];
            
            eval(strcat('neuron_animal.',NEURON{iNeuron},'.',strrep(ANIMAL{iAnimals},' ','_'),'.PMID = ','PMIDs',';'));
            
            clear PMIDs tmp_PMID
            
        end
        
    end

end

% Neuron & Animal & Phenomenon
for iNeuron=1:nNeuron
    
    for iAnimals=1:nAnimals
        
        for iPhenomenon=1:nPhenomenon
            
            disp(strcat('query-lucene-neuron-animals-phenomenon-',strrep(ANIMAL{iAnimals},' ','_'),'-',NEURON{iNeuron},'-',PHENOMENON{iPhenomenon},'.txt'));

            %%% data = loadjson(strcat('query-lucene-animals-neuron-phenomenon-',ANIMAL{iAnimals},'-',NEURON{iNeuron},'-',PHENOMENON{iPhenomenon},'.txt'));
            
            idx = strmatch(strcat('query-lucene-animals-neuron-phenomenon-',ANIMAL{iAnimals},'-',NEURON{iNeuron},'-',PHENOMENON{iPhenomenon}),labels,'exact');
            %%% Index = find(not(cellfun('isempty', idx)));
        
            nCount = all_queries(idx).result.response.numFound;
            
            if nCount > 0
        
                eval(strcat('neuron_animal_phenomenon.',NEURON{iNeuron},'.',strrep(ANIMAL{iAnimals},' ','_'),'.',PHENOMENON{iPhenomenon},'.count = ',int2str(nCount),';'));
                
                eval(strcat('all.neuron_animal_phenomenon.',NEURON{iNeuron},'.',strrep(ANIMAL{iAnimals},' ','_'),'.',PHENOMENON{iPhenomenon},'=',int2str(nCount),';'));
                
                tmp_PMID = all_queries(idx).result.facet_counts.facet_fields.PMID;
        
                nPMIDs = length(tmp_PMID);
                count = 0;
                for iPMID=1:2:nPMIDs
                    count = count + 1;
                    PMIDs{count} = tmp_PMID{iPMID};
                end
                
                allPMID = [allPMID,PMIDs];
                
                eval(strcat('neuron_animal_phenomenon.',NEURON{iNeuron},'.',strrep(ANIMAL{iAnimals},' ','_'),'.',PHENOMENON{iPhenomenon},'.PMID = ','PMIDs',';'));
                
                clear PMIDs tmp_PMID
                
            end
            
        end
        
    end

end

%%% SEE WHICH PMIDs ARE OUTSIDE ALL QUERIES

%%% data = loadjson('query-lucene-Retina-healthy-PMIDs.txt');

idx = strmatch('query-lucene-Retina-healthy-PMIDs',labels,'exact');
%%% Index = find(not(cellfun('isempty', idx)));
            
tmp_PMID = all_queries(idx).result.facet_counts.facet_fields.PMID;

nPMIDs = length(tmp_PMID);
count = 0;
for iPMID=1:2:nPMIDs
    count = count + 1;
    PMIDs{count} = tmp_PMID{iPMID};
end


allPMID = unique(allPMID);

idx = ismember(PMIDs,allPMID);
idx = ~idx;
idx = find(idx);
unknown = PMIDs(idx); %%% UNKNOWN PMIDs

save(strcat('__PMIDs-',goal,'_v',version,'.mat'),'all','animal','neuron','phenomenon','animal_neuron','animal_neuron_phenomenon','neuron_animal','neuron_animal_phenomenon','unknown','allPMID');

eval(strcat(goal,'.all = ','all;'));

eval(strcat('dataxml = struct2xml(',goal,');'));

fid = fopen(strcat('__PMIDs-',goal,'_v',version,'.xml'),'w');
fprintf(fid,'%s\n',dataxml);
fclose(fid);

