function processBigXMLv4(xmlfile)

%xmlfile = 'attention.xml';

%%% READS XML-PUBMED DATA FILE AND EXTRACT KEYWORDS DATA TO DO STATISTICS

xml = fopen(xmlfile, 'r', 'l', 'UTF-8');

expectedLine{1} = '<?xml version="1.0"?>';
    
expectedLine{2} = '<PubmedArticleSet>';

expectedLine{3} = '</PubmedArticleSet>';
    
expectedLine{4} = '<PubmedArticle>';

expectedLine{5} = '</PubmedArticle>';
    
expectedLine{6} = '<KeywordList';

expectedLine{7} = '</KeywordList>';
     
expectedLine{8} = '<Keyword';

expectedLine{9} = '</Keyword>';
      
expectedLine{10} = '<Journal>';

expectedLine{11} = '</Journal>';

expectedLine{12} = '<Title>'; %% Journal Title
      
expectedLine{13} = '</Title>';

expectedLine{14} = '<Abstract>';

expectedLine{15} = '</Abstract>';

expectedLine{16} = '<AuthorList>';

expectedLine{17} = '</AuthorList>';

expectedLine{18} = '<Author>';

expectedLine{19} = '</Author>';

expectedLine{20} = '<ArticleTitle>';

expectedLine{21} = '</ArticleTitle>';

PubMed = struct('ArticleTitle',char.empty,'Abstract',char.empty,'JournalTitle',char.empty,'Authors',char.empty,'Keywords',char.empty);

wordNet = load('adj.mat');
adj = strtrim(wordNet.uniqueWords);

wordNet = load('adv.mat');
adv = strtrim(wordNet.uniqueWords);

wordNet = load('verb.mat');
verb = strtrim(wordNet.uniqueWords);

wordNet = load('noun-pl.mat');
noun = strtrim(wordNet.allNouns);

someWords = load('someToRemove.mat');
%someWordsToRemove = strtrim(someWords.someToRemove);
someWordsToRemove = char.empty; %% just to test without removing words beyond wordnet

nArticle = 0;
nJournalTitle = 0;
nKeywords = 0;
nAuthors = 0;
nArticleTitle = 0;
nAbstract = 0;

l = 0;
step = 100000;
threshold = 100000;

while ~feof(xml)
    
    line = fgetl(xml);
    l = l + 1;
    
    if l > threshold; 
        
        disp(strcat('start loop line:',int2str(l))); 
        threshold = threshold + step;
    
    end
        
    lineToBeFound = '<PubmedArticle>';
    idx_start_of_article_string = strfind(line,lineToBeFound);
    
    lineToBeFound = '<PubmedBookArticle>';
    idx_start_of_book_article_string = strfind(line,lineToBeFound);
    
    if ~isempty(idx_start_of_article_string) || ~isempty(idx_start_of_book_article_string)
        
        nArticle = nArticle + 1;
            
        isNotEndOfArticle = true;
        
        %if ~isempty(idx_start_of_article_string); disp('<PubmedArticle>'); end
        %if ~isempty(idx_start_of_book_article_string); disp('<PubmedBookArticle>'); end 
        
        while isNotEndOfArticle
            
            lineToBeFound = '<Title>'; %% Journal Title
            idx_title_closed_string = strfind(line,lineToBeFound);
            
            lineToBeFound = '<Title '; %% Journal Title
            idx_title_open_string = strfind(line,lineToBeFound);
    
            lineToBeFound = '<BookTitle>'; %% Journal Title
            idx_booktitle_closed_string = strfind(line,lineToBeFound);
    
            lineToBeFound = '<BookTitle '; %% Journal Title
            idx_booktitle_open_string = strfind(line,lineToBeFound);
        
            lineToBeFound = '<KeywordList'; %% Keywords
            idx_keywordlist_string = strfind(line,lineToBeFound);
             
            lineToBeFound = '<Abstract'; %% Abstract
            idx_abstract_string = strfind(line,lineToBeFound);
            
            lineToBeFound = '<AuthorList'; %% Authors
            idx_authorlist_string = strfind(line,lineToBeFound);
            
            lineToBeFound = '<Authors';
            idx_authors_string = strfind(line,lineToBeFound);
            
            lineToBeFound = '<ArticleTitle';
            idx_articletitle_string = strfind(line,lineToBeFound);
            
            lineToBeFound = '<CollectionTitle';
            idx_collectiontitle_string = strfind(line,lineToBeFound);
            
            lineToBeFound = '</PubmedArticle>';
            idx_end_of_article_string = strfind(line,lineToBeFound);
            
            lineToBeFound = '</PubmedBookArticle>';
            idx_end_of_book_article_string = strfind(line,lineToBeFound);
            
            lineToBeFound = '<AbstractText';
            idx_abstracttext_string = strfind(line,lineToBeFound);
           
            testEverything = [~isempty(idx_title_closed_string) ~isempty(idx_title_open_string) ~isempty(idx_booktitle_closed_string) ~isempty(idx_booktitle_open_string) ~isempty(idx_keywordlist_string) ~isempty(idx_abstract_string) ~isempty(idx_authorlist_string) ~isempty(idx_authors_string) ~isempty(idx_articletitle_string) ~isempty(idx_collectiontitle_string) ~isempty(idx_end_of_article_string) ~isempty(idx_end_of_book_article_string)];
            
            idx_string_found = find(testEverything);
            
            if ~isempty(idx_string_found)
            
                switch idx_string_found

                    case 1 %% Title (Journal)
                        
                        %disp('<Title>');
                        
                        nJournalTitle = nJournalTitle + 1;
                        
                        [PubMed(nArticle).JournalTitle, line] = getJournalTitle(xml,line,'<Title>','</Title>');
                        
                    case 2 %% Title (Journal)
                        
                        %disp('<Title ');
                        
                        nJournalTitle = nJournalTitle + 1;
                        
                        [PubMed(nArticle).JournalTitle, line] = getJournalTitle(xml,line,'<Title ','</Title>');
                        
                    case 3 %% Title (Journal)
                        
                        %disp('<BookTitle>');
                        
                        nJournalTitle = nJournalTitle + 1;
                        
                        [PubMed(nArticle).JournalTitle, line] = getJournalTitle(xml,line,'<BookTitle>','</BookTitle>');
                        
                    case 4 %% Title (Journal)
                        
                        %disp('<BookTitle ');
                        
                        nJournalTitle = nJournalTitle + 1;
                        
                        [PubMed(nArticle).JournalTitle, line] = getJournalTitle(xml,line,'<BookTitle ','</BookTitle>');

                    case 5 %% KeywordList
                        
                        %disp('<KeywordList');
                        
                        nKeywords = nKeywords + 1;

                        [PubMed(nArticle).Keywords, line] = getKeywords(xml,line);

                    case 6 %% Abstract
                        
                        if isempty(idx_abstracttext_string)
                            
                            %disp('<Abstract');

                            nAbstract = nAbstract + 1;

                            [PubMed(nArticle).Abstract, line] = getAbstract(xml,line,adj,adv,verb,noun,someWordsToRemove);
                            
                        else
                            
                            line = fgetl(xml);
                            l = l + 1;
                            
                        end

                    case 7 %% AuthorList
                        
                        %disp('<AuthorList');

                        nAuthors = nAuthors + 1;
        
                        [PubMed(nArticle).Authors, line] = getAuthors(xml,line,'</AuthorList>');
                        
                    case 8 %% Authors
                        
                        %disp('<Authors');

                        nAuthors = nAuthors + 1;
        
                        [PubMed(nArticle).Authors, line] = getAuthors(xml,line,'</Authors>');

                    case 9 %% ArticleTitle
                        
                        %disp('<ArticleTitle');
                        
                        nArticleTitle = nArticleTitle + 1;

                        [PubMed(nArticle).ArticleTitle, line] = getArticleTitle(xml,line,'<ArticleTitle','</ArticleTitle',adj,adv,verb,noun,someWordsToRemove);
                        
                    case 10 %% CollectionTitle
                        
                        %disp('<CollectionTitle');
                        
                        nArticleTitle = nArticleTitle + 1;

                        [PubMed(nArticle).ArticleTitle, line] = getArticleTitle(xml,line,'<CollectionTitle','<CollectionTitle',adj,adv,verb,noun,someWordsToRemove);

                    case 11 %% End Of Article
                        
                        %disp('</PubmedArticle>');

                        isNotEndOfArticle = false;
                        
                    case 12 %% End Of Article
                        
                        %disp('</PubmedBookArticle>');

                        isNotEndOfArticle = false;

                end
            
            end
            
            if isempty(idx_string_found)
                
                line = fgetl(xml);
                l = l + 1;
                
            end
            
        end
        
    end
    
end

statFileName = strcat(xmlfile(1:end-4),'-PubMed-v4.mat');

save(statFileName,'PubMed','nJournalTitle','nKeywords','nAbstract','nAuthors','nArticleTitle');

fclose('all');

return


function [title, fileline] = getJournalTitle(pubmedfile,fileline,start_section,end_section)
    
    title = [];
    
    idx_start_of_section_title = strfind(fileline,start_section);

    idx_end_of_section_title = strfind(fileline,end_section);

    if ~isempty(idx_start_of_section_title) && ~isempty(idx_end_of_section_title)

        start_ = idx_start_of_section_title+length(start_section);

        end_ = idx_end_of_section_title - 1;

        title = fileline(start_:end_);

    end    
    
    fileline = fgetl(pubmedfile);
    
return

function [keywords, fileline] = getKeywords(pubmedfile,fileline)

    FID = pubmedfile;
    POSITION = ftell(FID);
    
    idx_end_of_section_keywordlist = [];
    
    nKeys = 1;
    
    keywords = char.empty;
    
    while isempty(idx_end_of_section_keywordlist)
    
        end_of_section_keywordlist = '</KeywordList>';
        idx_end_of_section_keywordlist = strfind(fileline,end_of_section_keywordlist);
        
        start_of_keyword = '<Keyword ';
        idx_start_of_keyword = strfind(fileline,start_of_keyword);
        
        end_of_keyword = '</Keyword>';
        idx_end_of_keyword = strfind(fileline,end_of_keyword);
        
        if ~isempty(idx_start_of_keyword)
            
            end_key = '>';
            idx_end_key = strfind(fileline,end_key);
        
            keywords{nKeys} = fileline(idx_end_key(1)+1:idx_end_of_keyword-1);
        
            nKeys = nKeys + 1;
            
        end
        
        fileline = fgetl(pubmedfile);
        
    end
    
    if size(keywords,1) ~= 1; keywords = keywords'; end
    
    keywords = lower(keywords);
    
%     ORIGIN = 'bof';
%     OFFSET = POSITION + 1;
%     FID = pubmedfile;
%     STATUS = fseek(FID, OFFSET, ORIGIN);
     
return

function [abstract, fileline] = getAbstract(pubmedfile,fileline,adj,adv,verb,noun,someWordsToRemove)
        
    FID = pubmedfile;
    POSITION = ftell(FID);
    
    allAbstract = [];
    idx_start_of_section_abstract = [];
    idx_end_of_section_abstract = [];
    
    while isempty(idx_end_of_section_abstract)
    
        start_of_section_abstract = '<Abstract';
        idx_start_of_section_abstract = strfind(fileline,start_of_section_abstract);
    
        end_of_section_abstract = '</Abstract>';
        idx_end_of_section_abstract = strfind(fileline,end_of_section_abstract);
        
        start_of_sub_section_abstracttext = '<AbstractText';
        idx_start_of_sub_section_abstracttext = strfind(fileline,start_of_sub_section_abstracttext);
    
        if ~isempty(idx_start_of_sub_section_abstracttext)
    
            end_string = '>';
            idx_end = strfind(fileline,end_string);
            
            fileline = fileline(idx_end(1)+1:end);
            
        end
        
        if ~isempty(idx_start_of_section_abstract)
            
            start_ = idx_start_of_section_abstract+length(start_of_section_abstract)+1;
            
        else
            
            start_ = 1;
            
        end
        
        
        if ~isempty(idx_end_of_section_abstract); 
        
            end_ = idx_end_of_section_abstract - 1; 
    
        else
        
            end_ = length(fileline);
    
        end
            
        allAbstract = strcat(allAbstract,fileline(start_:end_));
        
        fileline = fgetl(pubmedfile);
        
    end

    allAbstract = removeUndesiredStrings(allAbstract);
    
    idx_spaces = find(isspace(allAbstract));
    
    if ~isempty(idx_spaces)
        
        abstractKeys{1} = allAbstract(1:idx_spaces(1)-1);
    
        for i=1:length(idx_spaces)-1
       
            abstractKeys{i+1} = allAbstract(idx_spaces(i)+1:idx_spaces(i+1)-1);
        
        end
    
        abstractKeys{length(idx_spaces)} = allAbstract(idx_spaces(end)+1:end);
        
        abstractKeys = strtrim(abstractKeys);
    
        [abstract, countAbstract] = count_unique(abstractKeys);
    
        abstract(strcmp('',abstract)) = [];
      
    else
        
        abstract = allAbstract;
        
    end
    
    nWordsInAbstract = length(abstract);
    
    idx_adj = [];
    idx_adv = [];
    idx_verb = [];
    idx_some = [];
    
    idx_words_to_remove = [];
    
    if isempty(idx_spaces)
        
        idx_adj = find(strmatch(lower(abstract(1:end-2)),lower(adj)));
        idx_adv = find(strmatch(lower(abstract(1:end-2)),lower(adv)));
        idx_verb = find(strmatch(lower(abstract(1:end-2)),lower(verb)));
        idx_some = find(strmatch(lower(abstract(1:end-2)),lower(someWordsToRemove)));
        
        if ~isempty(idx_adj) || ~isempty(idx_adv) || ~isempty(idx_verb) || ~isempty(idx_some)
            
            abstract = [];
            
        end
        
    else
    
        for iWords=1:nWordsInAbstract

            idx_adj = find(strmatch(lower(abstract{iWords}(1:end-2)),lower(adj)));
            idx_adv = find(strmatch(lower(abstract{iWords}(1:end-2)),lower(adv)));
            idx_verb = find(strmatch(lower(abstract{iWords}(1:end-2)),lower(verb)));
            idx_some = find(strmatch(lower(abstract{iWords}(1:end-2)),lower(someWordsToRemove)));

            if ~isempty(idx_adj) || ~isempty(idx_adv) || ~isempty(idx_verb) || ~isempty(idx_some)

                idx_words_to_remove = [idx_words_to_remove, iWords];

            end

        end

        if ~isempty(idx_words_to_remove)

            abstract(idx_words_to_remove) = [];

        end
        
    end
    
    if size(abstract,1) ~= 1; abstract = abstract'; end
    
    abstract = lower(abstract);
    
%     ORIGIN = 'bof';
%     OFFSET = POSITION + 1;
%     FID = pubmedfile;
%     STATUS = fseek(FID, OFFSET, ORIGIN);
    
return

function [authors, fileline] = getAuthors(pubmedfile,fileline,end_section) 
        
    FID = pubmedfile;
    POSITION = ftell(FID);
    
    idx_end_of_section_authorlist = [];
    
    iAuthors = 0;
    
    authors = cell.empty;
    
    while isempty(idx_end_of_section_authorlist)
        
        idx_end_of_section_authorlist = strfind(fileline,end_section);
        
        start_of_author_one = '<Author ';
        start_of_author_two = '<Author>';
        idx_start_of_author_one = strfind(fileline,start_of_author_one);
        idx_start_of_author_two = strfind(fileline,start_of_author_two);
        
        start_of_lastname = '<LastName>';
        idx_start_of_lastname = strfind(fileline,start_of_lastname);
        
        end_of_lastname = '</LastName>';
        idx_end_of_lastname = strfind(fileline,end_of_lastname);
        
        start_of_forename = '<ForeName>';
        idx_start_of_forename = strfind(fileline,start_of_forename);
        
        end_of_forename = '</ForeName>';
        idx_end_of_forename = strfind(fileline,end_of_forename);
        
        if ~isempty(idx_start_of_author_one) || ~isempty(idx_start_of_author_two)
            
            iAuthors = iAuthors + 1;
            
        end
        
        if ~isempty(idx_start_of_lastname)
        
            authors{iAuthors,1} = fileline(idx_start_of_lastname+length(start_of_lastname):idx_end_of_lastname-1);
            
        end
        
        if ~isempty(idx_start_of_forename)
        
            authors{iAuthors,2} = fileline(idx_start_of_forename+length(start_of_forename):idx_end_of_forename-1);
            
        end
        
        fileline = fgetl(pubmedfile);
        
    end
    
%     ORIGIN = 'bof';
%     OFFSET = POSITION + 1;
%     FID = pubmedfile;
%     STATUS = fseek(FID, OFFSET, ORIGIN);
    
return

function [articletitlekeywords, fileline] = getArticleTitle(pubmedfile,fileline,start_section,end_section,adj,adv,verb,noun,someWordsToRemove)
        
    FID = pubmedfile;
    POSITION = ftell(FID);
    
    articletitle = [];
    idx_start_of_section_title = [];
    idx_end_of_section_title = [];
    
    while isempty(idx_end_of_section_title)
    
        idx_start_of_section_title = strfind(fileline,start_section);
    
        idx_end_of_section_title = strfind(fileline,end_section);
    
        if ~isempty(idx_start_of_section_title)
           
            end_string = '>';
           
            idx_end = strfind(fileline,end_string);
            
            %fileline = fileline(idx_end(1)+1:end);
            
            %start_ = idx_start_of_section_title+length(start_of_section_title);
            
            start_ = idx_end(1)+1;
            
        else
            
            start_ = 1;
            
        end
        
        
        if ~isempty(idx_end_of_section_title); 
        
            end_ = idx_end_of_section_title - 1; 
    
        else
        
            end_ = length(fileline);
    
        end
    
        articletitle = strcat(articletitle,fileline(start_:end_));
        
        fileline = fgetl(pubmedfile);
        
    end
    
    articletitle = removeUndesiredStrings(articletitle);
   
    idx_spaces = find(isspace(articletitle));
    
    if ~isempty(idx_spaces)
        
        articletitleKeys{1} = articletitle(1:idx_spaces(1)-1);
    
        for i=1:length(idx_spaces)-1
       
            articletitleKeys{i+1} = articletitle(idx_spaces(i)+1:idx_spaces(i+1)-1);
        
        end
    
        articletitleKeys{length(idx_spaces)+1} = articletitle(idx_spaces(end)+1:end);
        
        articletitleKeys = strtrim(articletitleKeys);
    
        [articletitlekeywords, nArticletitlekeywords] = count_unique(articletitleKeys);

    else
        
        articletitlekeywords = articletitle;
        
    end
    
    nWordsInArticleTitle = length(articletitlekeywords);
    
    idx_adj = [];
    idx_adv = [];
    idx_verb = [];
    idx_some = [];
    
    idx_words_to_remove = [];
    
    if isempty(idx_spaces)
        
        idx_adj = find(strmatch(lower(articletitlekeywords(1:end-2)),lower(adj)));
        idx_adv = find(strmatch(lower(articletitlekeywords(1:end-2)),lower(adv)));
        idx_verb = find(strmatch(lower(articletitlekeywords(1:end-2)),lower(verb)));
        idx_some = find(strmatch(lower(articletitlekeywords(1:end-2)),lower(someWordsToRemove)));
        
        if ~isempty(idx_adj) || ~isempty(idx_adv) || ~isempty(idx_verb) || ~isempty(idx_some)
            
            articletitlekeywords = [];
            
        end
        
    else
        
        for iWords=1:nWordsInArticleTitle

            idx_adj = find(strmatch(lower(articletitlekeywords{iWords}(1:end-2)),lower(adj)));
            idx_adv = find(strmatch(lower(articletitlekeywords{iWords}(1:end-2)),lower(adv)));
            idx_verb = find(strmatch(lower(articletitlekeywords{iWords}(1:end-2)),lower(verb)));
            idx_some = find(strmatch(lower(articletitlekeywords{iWords}(1:end-2)),lower(someWordsToRemove)));

            if ~isempty(idx_adj) || ~isempty(idx_adv) || ~isempty(idx_verb) || ~isempty(idx_some)

                idx_words_to_remove = [idx_words_to_remove, iWords];

            end

        end

        if ~isempty(idx_words_to_remove)

            articletitlekeywords(idx_words_to_remove) = [];

        end
        
    end
    
    if size(articletitlekeywords,1) ~= 1; articletitlekeywords = articletitlekeywords'; end
    
    articletitlekeywords = lower(articletitlekeywords);
    
%     ORIGIN = 'bof';
%     OFFSET = POSITION + 1;
%     FID = pubmedfile;
%     STATUS = fseek(FID, OFFSET, ORIGIN);
        
return

function output_string = removeUndesiredStrings(input_string)
    
    input_string = strrep(input_string,'.','');
    input_string = strrep(input_string,',','');
    input_string = strrep(input_string,'; ','');
    input_string = strrep(input_string,'?','');
    input_string = strrep(input_string,'!','');
    input_string = strrep(input_string,'- ',' ');
    input_string = strrep(input_string,'`','');
    input_string = strrep(input_string,'''','');
    
    input_string = strrep(input_string,' in ',' ');
    input_string = strrep(input_string,' on ',' ');
    input_string = strrep(input_string,' at ',' ');
    
    input_string = strrep(input_string,'In ','');
    input_string = strrep(input_string,'On ','');
    input_string = strrep(input_string,'At ','');
    
    input_string = strrep(input_string,' for ',' ');
    input_string = strrep(input_string,' from ',' ');
    input_string = strrep(input_string,' to ',' ');
    input_string = strrep(input_string,' by ',' ');
    input_string = strrep(input_string,' of ',' ');
    input_string = strrep(input_string,' with ',' ');
    input_string = strrep(input_string,' within ',' ');
    input_string = strrep(input_string,' without ',' ');
    
    input_string = strrep(input_string,'For ','');
    input_string = strrep(input_string,'From ','');
    input_string = strrep(input_string,'To ','');
    input_string = strrep(input_string,'By ','');
    input_string = strrep(input_string,'Of ','');
    input_string = strrep(input_string,'With ','');
    input_string = strrep(input_string,'Within ','');
    input_string = strrep(input_string,'Without ','');
    
    input_string = strrep(input_string,' the ',' ');
    input_string = strrep(input_string,' a ',' ');
    input_string = strrep(input_string,' an ',' ');
    
    input_string = strrep(input_string,'The ','');
    input_string = strrep(input_string,'A ','');
    input_string = strrep(input_string,'An ','');
    
    input_string = strrep(input_string,' this ',' ');
    input_string = strrep(input_string,' that ',' ');
    input_string = strrep(input_string,' these ',' ');
    input_string = strrep(input_string,' those ',' ');
    
    input_string = strrep(input_string,'This ','');
    input_string = strrep(input_string,'That ','');
    input_string = strrep(input_string,'These ','');
    input_string = strrep(input_string,'Those ','');
    
    input_string = strrep(input_string,' before ',' ');
    input_string = strrep(input_string,' after ',' ');
    
    input_string = strrep(input_string,'Before ','');
    input_string = strrep(input_string,'After ','');
    
    input_string = strrep(input_string,' other ',' ');
    input_string = strrep(input_string,' another ',' ');
    
    input_string = strrep(input_string,'Other ','');
    input_string = strrep(input_string,'Another ','');
    
    input_string = strrep(input_string,' or ',' ');
    
    input_string = strrep(input_string,' then ',' ');
    input_string = strrep(input_string,' than ',' ');
    
    input_string = strrep(input_string,'Then ','');
    input_string = strrep(input_string,'Than ','');
    
    input_string = strrep(input_string,'I ','');
    input_string = strrep(input_string,'You ','');
    input_string = strrep(input_string,'He ','');
    input_string = strrep(input_string,'She ','');
    input_string = strrep(input_string,'It ','');
    input_string = strrep(input_string,'We ','');
    input_string = strrep(input_string,'They ','');
    
    input_string = strrep(input_string,' i ',' ');
    input_string = strrep(input_string,' you ',' ');
    input_string = strrep(input_string,' he ',' ');
    input_string = strrep(input_string,' she ',' ');
    input_string = strrep(input_string,' it ',' ');
    input_string = strrep(input_string,' we ',' ');
    input_string = strrep(input_string,' they ',' ');
    
    input_string = strrep(input_string,' my ',' ');
    input_string = strrep(input_string,' meine ',' ');
    input_string = strrep(input_string,' your ',' ');
    input_string = strrep(input_string,' yours ',' ');
    input_string = strrep(input_string,' his ',' ');
    input_string = strrep(input_string,' her ',' ');
    input_string = strrep(input_string,' its ',' ');
    input_string = strrep(input_string,' our ',' ');
    input_string = strrep(input_string,' ours ',' ');
    input_string = strrep(input_string,' their ',' ');
    
    input_string = strrep(input_string,'My ',' ');
    input_string = strrep(input_string,'Meine ',' ');
    input_string = strrep(input_string,'Your ',' ');
    input_string = strrep(input_string,'Yours ',' ');
    input_string = strrep(input_string,'His ',' ');
    input_string = strrep(input_string,'Her ',' ');
    input_string = strrep(input_string,'Its ',' ');
    input_string = strrep(input_string,'Our ',' ');
    input_string = strrep(input_string,'Ours ',' ');
    input_string = strrep(input_string,'Their ',' ');
    
    input_string = strrep(input_string,' me ',' ');
    input_string = strrep(input_string,' you ',' ');
    input_string = strrep(input_string,' him ',' ');
    %input_string = strrep(input_string,' her ',' ');
    %input_string = strrep(input_string,' it ',' ');
    input_string = strrep(input_string,' us ',' ');
    input_string = strrep(input_string,' them ',' ');
    
    input_string = strrep(input_string,'Me ',' ');
    input_string = strrep(input_string,'You ',' ');
    input_string = strrep(input_string,'Him ',' ');
    %input_string = strrep(input_string,'Her ',' ');
    %input_string = strrep(input_string,'It ',' ');
    %input_string = strrep(input_string,'Us ',' ');
    input_string = strrep(input_string,'Them ',' ');
    
    input_string = strrep(input_string,'Can ','');
    input_string = strrep(input_string,'Both ','');
    input_string = strrep(input_string,'But ','');
    
    input_string = strrep(input_string,' can ',' ');
    input_string = strrep(input_string,' both ',' ');
    input_string = strrep(input_string,' but ',' ');
    
    input_string = strrep(input_string,' only ',' ');
    input_string = strrep(input_string,' some ',' ');
    input_string = strrep(input_string,' most ',' ');
    input_string = strrep(input_string,' more ',' ');
    
    input_string = strrep(input_string,'Only ','');
    input_string = strrep(input_string,'Some ','');
    input_string = strrep(input_string,'Most ','');
    input_string = strrep(input_string,'More ','');
    
    input_string = strrep(input_string,' further ',' ');
    input_string = strrep(input_string,' farther ',' ');
    
    input_string = strrep(input_string,'Further ','');
    input_string = strrep(input_string,'Farther ','');
    
    input_string = strrep(input_string,' good ',' ');
    input_string = strrep(input_string,' well ',' ');
    input_string = strrep(input_string,' bad ',' ');
    
    input_string = strrep(input_string,'Good ',' ');
    input_string = strrep(input_string,'Well ',' ');
    input_string = strrep(input_string,'Bad ',' ');
    
    input_string = strrep(input_string,' which ',' ');
    input_string = strrep(input_string,' what ',' ');
    input_string = strrep(input_string,' why ',' ');
    input_string = strrep(input_string,' where ',' ');
    input_string = strrep(input_string,' when ',' ');
    input_string = strrep(input_string,' who ',' ');
    input_string = strrep(input_string,' while ',' ');
    
    input_string = strrep(input_string,'Which ','');
    input_string = strrep(input_string,'What ','');
    input_string = strrep(input_string,'Why ','');
    input_string = strrep(input_string,'Where ','');
    input_string = strrep(input_string,'When ','');
    input_string = strrep(input_string,'Who ','');
    input_string = strrep(input_string,'While ','');
    
    input_string = strrep(input_string,' start ',' ');
    input_string = strrep(input_string,' begin ',' ');
    input_string = strrep(input_string,' end ',' ');
    input_string = strrep(input_string,' finish ',' ');
    
    input_string = strrep(input_string,'Start ','');
    input_string = strrep(input_string,'Begin ','');
    input_string = strrep(input_string,'End ','');
    input_string = strrep(input_string,'Finish ','');
    
    input_string = strrep(input_string,'During ','');
    input_string = strrep(input_string,' during ',' ');
    
    input_string = strrep(input_string,'Until ','');
    input_string = strrep(input_string,' until ',' ');
    
    input_string = strrep(input_string,' above ',' ');
    input_string = strrep(input_string,' below ',' ');
    input_string = strrep(input_string,' beside ',' ');
    input_string = strrep(input_string,' out ',' ');
    input_string = strrep(input_string,' under ',' ');
    input_string = strrep(input_string,' between ',' ');
    input_string = strrep(input_string,' behind ',' ');
    input_string = strrep(input_string,' as ',' ');
    input_string = strrep(input_string,' outside ',' ');
    input_string = strrep(input_string,' inside ',' ');
    
    input_string = strrep(input_string,'Above ','');
    input_string = strrep(input_string,'Below ','');
    input_string = strrep(input_string,'Beside ','');
    input_string = strrep(input_string,'Out ','');
    input_string = strrep(input_string,'Under ','');
    input_string = strrep(input_string,'Between ','');
    input_string = strrep(input_string,'Behind ','');
    input_string = strrep(input_string,'As ','');
    input_string = strrep(input_string,'Outside ','');
    input_string = strrep(input_string,'Inside ','');
    
    input_string = strrep(input_string,' is ',' ');
    input_string = strrep(input_string,' are ',' ');
    input_string = strrep(input_string,' were ',' ');
    input_string = strrep(input_string,' was ',' ');
    
    input_string = strrep(input_string,' has ',' ');
    input_string = strrep(input_string,' have ',' ');
    input_string = strrep(input_string,' had ',' ');
    
    input_string = strrep(input_string,' do ',' ');
    input_string = strrep(input_string,' does ',' ');
    input_string = strrep(input_string,' did ',' ');
    input_string = strrep(input_string,' done ',' ');
    
    input_string = strrep(input_string,' uses ',' ');
    input_string = strrep(input_string,' using ',' ');
    
    input_string = strrep(input_string,'Uses ','');
    input_string = strrep(input_string,'Using ','');
    
    input_string = strrep(input_string,' also ',' ');
    input_string = strrep(input_string,' too ',' ');
    
    input_string = strrep(input_string,'Also ','');
    input_string = strrep(input_string,'Too ','');
    
    input_string = strrep(input_string,' and ',' ');
    input_string = strrep(input_string,'And ','');
    
    input_string = strrep(input_string,' all ',' ');
    input_string = strrep(input_string,'All ','');
    
    input_string = strrep(input_string,' everything ',' ');
    input_string = strrep(input_string,'Everything ','');
    
    input_string = strrep(input_string,') ',' ');
    input_string = strrep(input_string,'] ',' ');
    input_string = strrep(input_string,'} ',' ');
    
    input_string = strrep(input_string,'):','');
    input_string = strrep(input_string,':(','');
    
    input_string = strrep(input_string,' (',' ');
    input_string = strrep(input_string,' [',' ');
    input_string = strrep(input_string,' {',' ');
    
    input_string = strrep(input_string,'[','');
    input_string = strrep(input_string,']','');
    
    input_string = strrep(input_string,'[]','');
    
    input_string = strrep(input_string,'<CopyrightInformation>','');
    input_string = strrep(input_string,'</CopyrightInformation>','');
    input_string = strrep(input_string,'Copyright','');
    input_string = strrep(input_string,'</AbstractText>','');
    input_string = strrep(input_string,'</Information>','');
    
    vowels = {' a ', ' b ', ' c ', ' d ', ' e ', ' f ', ' g ', ' h ', ' i ', ' j ', ' k ', ' l ', ' m ', ' n ', ' o ', ' p ', ' q ', ' r ', ' s ', ' t ', ' u ', ' v ', ' x ', ' y ', ' z '};
    
    nVowels = length(vowels);
    
    for iVowel=1:nVowels
        
        input_string = strrep(lower(input_string),vowels{iVowel},' ');
        
    end
    
    output_string = input_string;
        
return


