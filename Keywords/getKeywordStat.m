function getKeywordStat(pubmedxmlfile)

pubmedxml = xmlread(pubmedxmlfile);

paperslist = pubmedxml.getElementsByTagName('PubmedArticle');

nPaperslist = paperslist.getLength;

keywordlist = pubmedxml.getElementsByTagName('KeywordList');

nKeywordlist = keywordlist.getLength;

allKeys = cell(size(1));

k = 1;
for n=1:nKeywordlist
   
    list = keywordlist.item(n-1).getElementsByTagName('Keyword');
    
    nKeyInList = list.getLength;
    
    for nn=1:nKeyInList
    
        word = list.item(nn-1).getFirstChild.getData;
        
        if n == 1 && nn == 1
            
            allKeys{1} = char(word);
            
        else
           
            allKeys = [allKeys char(word)];
            
        end
    
    end
    
end

[allKeys, numKeys] = count_unique(allKeys,'float');

nAllKeys = length(allKeys);

keysToGether = struct('keys',cell(1,nAllKeys));

for n=1:nKeywordlist
   
    word = cell(size(1));
    
    list = keywordlist.item(n-1).getElementsByTagName('Keyword');
    
    nKeyInList = list.getLength;
    
    for nk=1:nKeyInList
               
        word{nk} = char(list.item(nk-1).getFirstChild.getData);
              
    end
       
    for nnn=1:nAllKeys

       exist = strcmp(allKeys{nnn},word);

       if sum(exist) ~= 0

          this_list = word;

          this_list(find(exist)) = [];

          remainedKeys = this_list; 

          keysToGether(nnn).keys = [keysToGether(nnn).keys remainedKeys];

       end

    end
       
end

for n=1:nAllKeys
    
    if ~isempty(keysToGether(n).keys)
    
        [uniqueKeys, freqKeys] = count_unique(keysToGether(n).keys,'float');
    
        keysToGether(n).uniqueKeys = uniqueKeys;
    
        keysToGether(n).freqKeys = freqKeys;
        
    else
        
        keysToGether(n).uniqueKeys = cell(size(1));
    
        keysToGether(n).freqKeys = 0;
        
    end
        
end

keywordStat.allKeys = allKeys;
keywordStat.numKeys = numKeys;
keywordStat.keysToGether = keysToGether;
keywordStat.nAllPapers = nPaperslist;
keywordStat.nKeywordPapers = nKeywordlist;

save('keywordStat.mat','keywordStat');


%end

