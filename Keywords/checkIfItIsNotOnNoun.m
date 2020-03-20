
wordNet = load('noun-pl.mat');
noun = strtrim(wordNet.allNouns);

nWords = length(allArticleWords);
allIDX = [];

parfor iWord=1:nWords
    
    if mod(iWord,1000) == 0; disp(int2str(iWord/1000)); end
    
    idx = find(strcmpi(allArticleWords{iWord},noun));
    
    if ~isempty(idx); allIDX = [allIDX, idx]; end
    
end

notIncludedWords = noun(allIDX);