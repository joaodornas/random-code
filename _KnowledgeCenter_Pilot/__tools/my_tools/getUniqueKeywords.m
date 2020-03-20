
function [uNouns, uVerbs, uAdverbs, uAdjectives, missing ] = getUniqueKeywords(allKeywords)

allKeywords = strrep(allKeywords,'(',' ');
allKeywords = strrep(allKeywords,')',' ');
allKeywords = strrep(allKeywords,'"',' ');
allKeywords = strrep(allKeywords,':',' ');
allKeywords = strrep(allKeywords,',',' ');
allKeywords = strrep(allKeywords,']',' ');
allKeywords = strrep(allKeywords,'[',' ');
allKeywords = strrep(allKeywords,'{',' ');
allKeywords = strrep(allKeywords,'}',' ');
allKeywords = strrep(allKeywords,'...',' ');
allKeywords = strrep(allKeywords,'.',' ');
allKeywords = strrep(allKeywords,';',' ');
allKeywords = strrep(allKeywords,'--',' ');
allKeywords = strrep(allKeywords,'?',' ');
allKeywords = strrep(allKeywords,'!',' ');
allKeywords = strrep(allKeywords,'/',' ');
allKeywords = strrep(allKeywords,'\',' ');
allKeywords = strrep(allKeywords,'''',' ');
allKeywords = strrep(allKeywords,'~',' ');
allKeywords = strrep(allKeywords,'@',' ');
allKeywords = strrep(allKeywords,'#',' ');
allKeywords = strrep(allKeywords,'$',' ');
allKeywords = strrep(allKeywords,'%',' ');
allKeywords = strrep(allKeywords,'^',' ');
allKeywords = strrep(allKeywords,'&',' ');
allKeywords = strrep(allKeywords,'*',' ');
allKeywords = strrep(allKeywords,'+',' ');
allKeywords = strrep(allKeywords,'|',' ');
allKeywords = strrep(allKeywords,'>',' ');
allKeywords = strrep(allKeywords,'<',' ');
allKeywords = strrep(allKeywords,'1',' ');
allKeywords = strrep(allKeywords,'2',' ');
allKeywords = strrep(allKeywords,'3',' ');
allKeywords = strrep(allKeywords,'4',' ');
allKeywords = strrep(allKeywords,'5',' ');
allKeywords = strrep(allKeywords,'6',' ');
allKeywords = strrep(allKeywords,'7',' ');
allKeywords = strrep(allKeywords,'8',' ');
allKeywords = strrep(allKeywords,'9',' ');
allKeywords = strrep(allKeywords,'0',' ');

nKeywords = length(allKeywords);

disp('splitting...');

all_these_keys = [];
for iKey=1:nKeywords
    
    this_key = allKeywords{iKey};
    these_keys = strsplit(this_key);
    these_keys = strtrim(these_keys);
    these_keys(strcmp('',these_keys)) = [];
    
    all_these_keys = [all_these_keys,these_keys];
    
end

allKeywords = lower(all_these_keys);
allKeywords = unique(allKeywords);

disp('getting adjectives...');

wordNet = load('adj.mat');
idx_found_adj = seeIfIKnowThisKeyword(allKeywords,wordNet);
if ~isempty(idx_found_adj)
    uAdjectives = allKeywords(idx_found_adj);
else
    uAdjectives = cell.empty;
end

disp('getting adverbs...');

wordNet = load('adv.mat');
idx_found_adv = seeIfIKnowThisKeyword(allKeywords,wordNet);
if ~isempty(idx_found_adv)
    uAdverbs = allKeywords(idx_found_adv);
else
    uAdverbs = cell.empty;
end

disp('getting verbs...');

wordNet = load('verb.mat');
idx_found_verb = seeIfIKnowThisKeyword(allKeywords,wordNet);
if ~isempty(idx_found_verb)
    uVerbs = allKeywords(idx_found_verb);
else
    uVerbs = cell.empty;
end

disp('getting nouns...');

wordNet = load('noun.mat');
idx_found_noun = seeIfIKnowThisKeyword(allKeywords,wordNet);
if ~isempty(idx_found_noun)
    uNouns = allKeywords(idx_found_noun);
else
    uNouns = cell.empty;
end

disp('getting missing ones...');

missing = allKeywords;
missing([idx_found_adj,idx_found_adv,idx_found_verb,idx_found_noun]) = [];

end

function idx_found = seeIfIKnowThisKeyword(allKeywords,wordNet)

idx_found = [];

wordNet = lower(strtrim(wordNet.uniqueWords));

for iKey=1:length(allKeywords)
    
    idx = find(strcmpi(allKeywords{iKey},wordNet));
    
    if isempty(idx)
        
        idx_found(iKey) = 0;
        
    else
        
        idx_found(iKey) = 1;
        
    end
    
end

idx_found = find(idx_found);

end
