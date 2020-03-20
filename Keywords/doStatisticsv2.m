tic

%load('attention-PubMed.mat');

%articleWords = [PubMed(1:end).ArticleTitle, PubMed(1:end).Abstract, PubMed(1:end).Keywords];

nArticles = length(PubMed);

articleWords = [];

for i=1:nArticles
    
    articleWords = [articleWords, PubMed(i).allKeywords'];
    
end
    
[allArticleWords, countAllArticleWords] = count_unique(articleWords);

[x,i] = sort(countAllArticleWords,'descend');

allArticleWords = allArticleWords(i);
countAllArticleWords = countAllArticleWords(i);

allArticleWordsSorted = sort(allArticleWords);

toc
    
    
    