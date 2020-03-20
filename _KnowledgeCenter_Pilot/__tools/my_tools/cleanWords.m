function words_cleaned = cleanWords( words, useWordNetsoft, useWordNethard, useDescriptive, descriptiveLabel, descriptiveVersion )

UICODE = 23;
minLength = 4;

if useWordNetsoft
    
    disp('getting adjectives...');
    adj = load('adj.mat');   

    disp('getting adverbs...');
    adv = load('adv.mat');

end

if useWordNethard
    
    disp('getting nouns...');
    noun = load('noun.mat');   

    disp('getting verbs...');
    verb = load('verb.mat');

end

if useDescriptive
   
    load(strcat('__descriptives-',descriptiveLabel,'_v',int2str(descriptiveVersion)));
    
end

nWords = length(words);
iiWord = 0;

for iWord=1:nWords
    
    word = words{iWord};
    
    word = lower(word);

    word = strrep(word,'(',' ');
    word = strrep(word,')',' ');
    word = strrep(word,'"',' ');
    word = strrep(word,':',' ');
    word = strrep(word,',',' ');
    word = strrep(word,']',' ');
    word = strrep(word,'[',' ');
    word = strrep(word,'{',' ');
    word = strrep(word,'}',' ');
    word = strrep(word,'...',' ');
    word = strrep(word,'.',' ');
    word = strrep(word,';',' ');
    word = strrep(word,'--',' ');
    word = strrep(word,'?',' ');
    word = strrep(word,'!',' ');
    word = strrep(word,'/',' ');
    word = strrep(word,'\',' ');
    word = strrep(word,'''',' ');
    word = strrep(word,'~',' ');
    word = strrep(word,'@',' ');
    word = strrep(word,'#',' ');
    word = strrep(word,'$',' ');
    word = strrep(word,'%',' ');
    word = strrep(word,'^',' ');
    word = strrep(word,'&',' ');
    word = strrep(word,'*',' ');
    word = strrep(word,'+',' ');
    word = strrep(word,'=',' ');
    word = strrep(word,'-',' ');
    word = strrep(word,'|',' ');
    word = strrep(word,'>',' ');
    word = strrep(word,'<',' ');
    word = strrep(word,'1',' ');
    word = strrep(word,'2',' ');
    word = strrep(word,'3',' ');
    word = strrep(word,'4',' ');
    word = strrep(word,'5',' ');
    word = strrep(word,'6',' ');
    word = strrep(word,'7',' ');
    word = strrep(word,'8',' ');
    word = strrep(word,'9',' ');
    word = strrep(word,'0',' ');
    word = strrep(word,'©',' ');
    word = strrep(word,'¢',' ');
    word = strrep(word,'°',' ');
    word = strrep(word,'â',' ');
    word = strrep(word,'á',' ');
    word = strrep(word,'à',' ');
    word = strrep(word,'ã',' ');
    word = strrep(word,'ô',' ');
    word = strrep(word,'ó',' ');
    word = strrep(word,'ò',' ');
    word = strrep(word,'õ',' ');
    word = strrep(word,'ö',' ');
    word = strrep(word,'ê',' ');
    word = strrep(word,'é',' ');
    word = strrep(word,'è',' ');
    word = strrep(word,'í',' ');
    word = strrep(word,'ì',' ');
    word = strrep(word,'ð',' ');
    word = strrep(word,'ÿ',' ');
    word = strrep(word,'€',' ');
    word = strrep(word,'¢',' ');
    word = strrep(word,'¡',' ');
    word = strrep(word,'¤',' ');
    word = strrep(word,'£',' ');
    word = strrep(word,'²',' ');
    word = strrep(word,'³',' ');
    word = strrep(word,'µ',' ');
    word = strrep(word,'¼',' ');
    word = strrep(word,'î',' ');
    word = strrep(word,'±',' ');
    
    word = strrep(word,' a ',' ');
    word = strrep(word,' an ',' ');
    word = strrep(word,' and ',' ');
    word = strrep(word,' are ',' ');
    word = strrep(word,' as ',' ');
    word = strrep(word,' at ',' ');
    word = strrep(word,' be ',' ');
    word = strrep(word,' but ',' ');
    word = strrep(word,' by ',' ');
    word = strrep(word,' for ',' ');
    word = strrep(word,' if ',' ');
    word = strrep(word,' in ',' ');
    word = strrep(word,' into ',' ');
    word = strrep(word,' is ',' ');
    word = strrep(word,' it ',' ');
    word = strrep(word,' no ',' ');
    word = strrep(word,' not ',' ');
    word = strrep(word,' of ',' ');
    word = strrep(word,' on ',' ');
    word = strrep(word,' or ',' ');
    word = strrep(word,' such ',' ');
    word = strrep(word,' that ',' ');
    word = strrep(word,' the ',' ');
    word = strrep(word,' their ',' ');
    word = strrep(word,' then ',' ');
    word = strrep(word,' there ',' ');
    word = strrep(word,' these ',' ');
    word = strrep(word,' they ',' ');
    word = strrep(word,' this ',' ');
    word = strrep(word,' to ',' ');
    word = strrep(word,' was ',' ');
    word = strrep(word,' will ',' ');
    word = strrep(word,' with ',' ');
    
    word = strrep(word,' a ',' ');
    word = strrep(word,' b ',' ');
    word = strrep(word,' c ',' ');
    word = strrep(word,' d ',' ');
    word = strrep(word,' e ',' ');
    word = strrep(word,' f ',' ');
    word = strrep(word,' g ',' ');
    word = strrep(word,' h ',' ');
    word = strrep(word,' i ',' ');
    word = strrep(word,' j ',' ');
    word = strrep(word,' k ',' ');
    word = strrep(word,' l ',' ');
    word = strrep(word,' m ',' ');
    word = strrep(word,' n ',' ');
    word = strrep(word,' o ',' ');
    word = strrep(word,' p ',' ');
    word = strrep(word,' r ',' ');
    word = strrep(word,' s ',' ');
    word = strrep(word,' t ',' ');
    word = strrep(word,' u ',' ');
    word = strrep(word,' v ',' ');
    word = strrep(word,' x ',' ');
    word = strrep(word,' y ',' ');
    word = strrep(word,' z ',' ');
    
    
    tmp_words = strsplit(word);
    tmp_words = strtrim(tmp_words);
    tmp_words(strcmp('',tmp_words)) = [];
    
    nTmp = length(tmp_words);
    
    for iTmp=1:nTmp
    
        tmp_word = tmp_words{iTmp};
        
        tmp_word = unicode2native(tmp_word);
        tmp_word(tmp_word<UICODE) = [];
        tmp_word = native2unicode(tmp_word);
        
        tmp_word = strtrim(tmp_word);
        tmp_word(strcmp(' ',tmp_word)) = [];
        
        if useWordNetsoft
        
            if length(tmp_word) >= minLength & isempty(find(strcmp(tmp_word,adj.uniqueWords))) & isempty(find(strcmp(tmp_word,adv.uniqueWords)))

                iiWord = iiWord + 1;

                words_cleaned{iiWord} = tmp_word;

            end
        
        elseif useWordNethard
        
            if length(tmp_word) >= minLength & isempty(find(strcmp(tmp_word,adj.uniqueWords))) & isempty(find(strcmp(tmp_word,adv.uniqueWords))) & isempty(find(strcmp(tmp_word,noun.uniqueWords))) & isempty(find(strcmp(tmp_word,verb.uniqueWords)))

                iiWord = iiWord + 1;

                words_cleaned{iiWord} = tmp_word;

            end
        
        elseif useWordNethard & useDescriptive
            
            if length(tmp_word) >= minLength & isempty(find(strcmp(tmp_word,adj.uniqueWords))) & isempty(find(strcmp(tmp_word,adv.uniqueWords))) & isempty(find(strcmp(tmp_word,noun.uniqueWords))) & isempty(find(strcmp(tmp_word,verb.uniqueWords))) & isempty(find(strcmp(tmp_word,animals))) & isempty(find(strcmp(tmp_word,neurons))) & isempty(find(strcmp(tmp_word,phenomenon)))

                iiWord = iiWord + 1;

                words_cleaned{iiWord} = tmp_word;

            end
            
        end
        
    end

end

end

