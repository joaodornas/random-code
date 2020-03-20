
tic

wordNet = load('noun.mat');
noun = lower(strtrim(wordNet.uniqueWords));

idx_not_found = [];
for i=1:length(allArticleWords)
    
    idx = [];
    idx = find(strcmpi(allArticleWords{i},noun));
    
    if isempty(idx)
        
        idx_not_found = [idx_not_found, i];
        
    end
    
end

words_not_noun = allArticleWords(idx_not_found);
    
toc
    
    