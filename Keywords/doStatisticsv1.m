
%load('attention-PubMed.mat');

nArticles = length(PubMed);
nPossibleWords = nArticles * 20;

allWords = char.empty;
countAllWords = nan(1,nPossibleWords);
journals = struct('names',char.empty,'count',[]);

for iArticle=1:nArticles
    
    if iArticle == 1; tic; end
    
    if mod(iArticle,1000) == 0; disp(int2str(iArticle)); toc; tic; end
    
    articleWords = char.empty;
    
    idx = [];
   
    idx = find([~isempty(PubMed(iArticle).ArticleTitle), ~isempty(PubMed(iArticle).Abstract), ~isempty(PubMed(iArticle).Keywords)]);
    
    if length(idx) == 1
    
        switch idx
            
            case 1
                
                articleWords = PubMed(iArticle).ArticleTitle;
                
            case 2
                
                articleWords = PubMed(iArticle).Abstract;
                
            case 3
                
                articleWords = PubMed(iArticle).Keywords;
                
        end
        
    else
        
        articleWords = [PubMed(iArticle).ArticleTitle, PubMed(iArticle).Abstract, PubMed(iArticle).Keywords];

    end
    
    if isempty(articleWords)
        
        uniArticleWords = [];
        
    else
        
        if ~iscellstr(articleWords); articleWords = cellstr(articleWords); end
    
        if length(articleWords) > 1; [uniArticleWords nuniArticleWords] = count_unique(articleWords); else uniArticleWords = char.empty; end
    
    end
    
    for iWords=1:length(uniArticleWords)
        
        word = uniArticleWords{iWords};
        
        idx_word_string = find(strcmpi(allWords,word));
        
        if isempty(idx_word_string)
            
            idx_now = length(allWords) + 1;
            
            allWords{idx_now} = word;
           
            countAllWords(idx_now) = 1;
            
            idx_word_string = idx_now;
            
        else
            
            countAllWords(idx_word_string) = countAllWords(idx_word_string) + 1;
            
        end
        
        journal = PubMed(iArticle).JournalTitle;
        
        if length(journals) >= idx_word_string

            idx_journal_string = find(strcmpi(journals(idx_word_string).names,journal));

            if isempty(idx_journal_string)

                idx_now = length(journals(idx_word_string).names) + 1;

                journals(idx_word_string).names{idx_now} = journal;

                journals(idx_word_string).count(idx_now) = 1;

                idx_journal_string = idx_now;

            else

                journals(idx_word_string).count(idx_journal_string) = journals(idx_word_string).count(idx_journal_string) + 1;

            end
            
        else
            
            journals(idx_word_string).names{1} = journal;
            
            journals(idx_word_string).count(1) = 1;
            
        end
        
    end
    
    clear articleWords
    clear uniArticleWords
    clear idx_word_string
    clear idx_journal_string
    clear idx_now
    clear word
        
end

toc




    
    