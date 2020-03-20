function [WordsUnique, nWordsUnique, frequencies, CondWordsUnique, nCondWordsUnique, CondFrequencies] = gpuGetWordFrequencies(WordMatrix)
   
nTotalWords = size(WordMatrix,1) * size(WordMatrix,2);

ConditionalWords = size(WordMatrix,2);

all_words = reshape(WordMatrix,1,[]);

[WordsUnique, nWordsUnique] = count_unique(all_words);

frequencies = nWordsUnique ./ nTotalWords;

for i=1:size(WordMatrix,2)
    
    bin_group = WordMatrix(:,i);
    
    [CondWordsUnique(i).words, nCondWordsUnique(i).nWords] = count_unique(bin_group);
    
    CondFrequencies(i).freq = nCondWordsUnique(i).nWords./ size(WordMatrix,1) ;
    
end

end

