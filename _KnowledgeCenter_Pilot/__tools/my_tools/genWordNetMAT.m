

folder = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_CODES/Visual-Attention-Awareness/my-toolbox/_TOOLBOX/WordNet-3.0/dict';

wordNetFileExc{1} = 'adj.exc';
wordNetFileExc{2} = 'adv.exc';
wordNetFileExc{3} = 'noun.exc';
wordNetFileExc{4} = 'verb.exc';

wordNetFileIndex{1} = 'index.adj';
wordNetFileIndex{2} = 'index.adv';
wordNetFileIndex{3} = 'index.noun';
wordNetFileIndex{4} = 'index.verb';

nWordNetFile = length(wordNetFileExc);

for iWordNet=1:nWordNetFile
    
   datafile = fopen(strcat(folder,'/',wordNetFileExc{iWordNet}), 'r', 'l', 'UTF-8');
   
   l = 0;
   
   allWords = char.empty;
   
   while ~feof(datafile)
       
        line = fgetl(datafile);
        l = l + 1;
        
        idx_spaces = find(isspace(line));
        
        if ~isempty(idx_spaces)
        
            words{1} = line(1:idx_spaces(1)-1);
    
            if length(idx_spaces) == 1
                
                words{2} = line(idx_spaces(1)+1:end);
                
            elseif length(idx_spaces) > 1
                
                for i=1:length(idx_spaces)-1
       
                    words{i+1} = line(idx_spaces(i)+1:idx_spaces(i+1)-1);
        
                end
            
                words{length(idx_spaces)} = line(idx_spaces(end)+1:end);
    
            end
    
        else
        
            words{1} = line;
        
        end
    
        allWords = [allWords, words];
        
        clear words;
        
   end
   
   fclose('all');
   
   datafile = fopen(strcat(folder,'/',wordNetFileIndex{iWordNet}), 'r', 'l', 'UTF-8');
   
   l = 0;
   
   while ~feof(datafile)
       
        line = fgetl(datafile);
        l = l + 1;
        
        idx_spaces = find(isspace(line));
        
        words = line(1:idx_spaces(1)-1);
        
        allWords = [allWords, words];
        
        clear words;
        
   end
   
   fclose('all');
   
   [uniqueWords, countWords] = count_unique(allWords);
    
   uniqueWords(strcmp('',uniqueWords)) = [];
   
   uniqueWords = unique(uniqueWords);
    
   save(strcat(wordNetFileExc{iWordNet}(1:end-4),'.mat'),'uniqueWords');
      
end