
prefix = 'retina';

animals = {'monkey', 'macaque', 'human', 'mammalian'};

keywords = {'rod', 'cone', 'photoreceptor', 'horizontal', 'bipolar', 'amacrine', 'ganglion', 'adaptation'};

nAnimals = length(animals);
nKeywords = length(keywords);

for iAnimal=1:nAnimals
    
    for iKeyword=1:nKeywords
        
        search_string = strcat(prefix,'+',animals{iAnimal},'+',keywords{iKeyword});
        
        %status = unix(sprintf('perl %s',strcat(search_string,'.pl')));
        
        status = unix(sprintf('med2xml %s.xml > %s.mod',search_string,search_string));
        
        status = unix(sprintf('xml2bib %s.mod > %s.bib',search_string,search_string));
        
    end
    
end