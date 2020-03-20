function processBigXMLv9(xmlfile)

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

PubMed = struct('allKeywords',char.empty,'JournalTitle',char.empty,'Authors',char.empty);

wordNet = load('adj.mat');
adj = lower(strtrim(wordNet.uniqueWords));

wordNet = load('adv.mat');
adv = lower(strtrim(wordNet.uniqueWords));

wordNet = load('verb.mat');
verb = lower(strtrim(wordNet.uniqueWords));

% wordNet = load('noun.mat');
% noun = lower(strtrim(wordNet.uniqueWords));
noun = char.empty;

% wordNet = load('noun-pl-sing.mat');
% plnoun = strtrim(wordNet.plNouns);
plnoun = char.empty;

nArticle = 0;
nJournalTitle = 0;
nKeywords = 0;
nAuthors = 0;
nArticleTitle = 0;
nAbstract = 0;

ArticleTitle = char.empty;
Abstract = char.empty;
Keywords = char.empty;

foundSection = false;

idx_start_of_article_string = [];
idx_start_of_book_article_string = [];
idx_title_closed_string = [];
idx_title_open_string = [];
idx_booktitle_closed_string = [];
idx_booktitle_open_string = [];
idx_keywordlist_string = [];
idx_abstract_string = [];
idx_authorlist_string = [];
idx_authors_string = [];
idx_articletitle_string = [];
idx_collectiontitle_string = [];
idx_end_of_article_string = [];
idx_end_of_book_article_string = [];
idx_abstracttext_string = [];

l = 0;
step = 100000;
threshold = 100000;

start = GetSecs;

% tic

while ~feof(xml)
    
    line = fgetl(xml);
    l = l + 1;
    foundSection = false;
    
    if l > threshold; 
        
        disp(strcat('start loop line:',int2str(l))); 
        threshold = threshold + step;
    
    end
        
    lineToBeFound = '<PubmedArticle>';
    if ~foundSection; idx_start_of_article_string = strfind(line,lineToBeFound); end
    if ~isempty(idx_start_of_article_string); foundSection = true; tic; end
    
    lineToBeFound = '<PubmedBookArticle>';
    if ~foundSection; idx_start_of_book_article_string = strfind(line,lineToBeFound); end
    if ~isempty(idx_start_of_book_article_string); foundSection = true; tic; end
    
    if ~isempty(idx_start_of_article_string) || ~isempty(idx_start_of_book_article_string)
        
        nArticle = nArticle + 1;
        
        if mod(nArticle,1000) == 0 && nArticle > 1
            
            nArticle
%             toc
%             tic
            
        end
        
        if nArticle > 1
            
            allKeywords = [ArticleTitle, Abstract, processedKeywords];
      
            if ~isempty(allKeywords)
            
                [allKeywords, nallKeywords] = count_unique(allKeywords);
                
            end
            
            PubMed(nArticle - 1).allKeywords = allKeywords;
            PubMed(nArticle - 1).keywords = Keywords;
            
            allKeywords = char.empty;
            ArticleTitle = char.empty;
            Abstract = char.empty;
            Keywords = char.empty;
            processedKeywords = char.empty;
           
            %break;
            
        end
            
        isNotEndOfArticle = true;
        
        while isNotEndOfArticle
            
            idx_start_of_article_string = [];
            idx_start_of_book_article_string = [];
            idx_title_closed_string = [];
            idx_title_open_string = [];
            idx_booktitle_closed_string = [];
            idx_booktitle_open_string = [];
            idx_keywordlist_string = [];
            idx_abstract_string = [];
            idx_authorlist_string = [];
            idx_authors_string = [];
            idx_articletitle_string = [];
            idx_collectiontitle_string = [];
            idx_end_of_article_string = [];
            idx_end_of_book_article_string = [];
            idx_abstracttext_string = [];
            
            lineToBeFound = '<Title>'; %% Journal Title
            if ~foundSection; idx_title_closed_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_title_closed_string); foundSection = true; end
            
            lineToBeFound = '<Title '; %% Journal Title
            if ~foundSection; idx_title_open_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_title_open_string); foundSection = true; end
                
            lineToBeFound = '<BookTitle>'; %% Journal Title
            if ~foundSection; idx_booktitle_closed_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_booktitle_closed_string); foundSection = true; end
                
            lineToBeFound = '<BookTitle '; %% Journal Title
            if ~foundSection; idx_booktitle_open_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_booktitle_open_string); foundSection = true; end
                
            lineToBeFound = '<KeywordList'; %% Keywords
            if ~foundSection; idx_keywordlist_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_keywordlist_string); foundSection = true; end
                
            lineToBeFound = '<Abstract'; %% Abstract
            if ~foundSection; idx_abstract_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_abstract_string); foundSection = true; end
                
            lineToBeFound = '<AuthorList'; %% Authors
            if ~foundSection; idx_authorlist_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_authorlist_string); foundSection = true; end
                
            lineToBeFound = '<Authors';
            if ~foundSection; idx_authors_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_authors_string); foundSection = true; end
                
            lineToBeFound = '<ArticleTitle';
            if ~foundSection; idx_articletitle_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_articletitle_string); foundSection = true; end
                
            lineToBeFound = '<CollectionTitle';
            if ~foundSection; idx_collectiontitle_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_collectiontitle_string); foundSection = true; end
                
            lineToBeFound = '</PubmedArticle>';
            if ~foundSection; idx_end_of_article_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_end_of_article_string); foundSection = true; end
                
            lineToBeFound = '</PubmedBookArticle>';
            if ~foundSection; idx_end_of_book_article_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_end_of_book_article_string); foundSection = true; end
                
            lineToBeFound = '<AbstractText';
            if ~foundSection; idx_abstracttext_string = strfind(line,lineToBeFound); end
            if ~isempty(idx_abstracttext_string); foundSection = true; end
                
            testEverything = [~isempty(idx_title_closed_string) ~isempty(idx_title_open_string) ~isempty(idx_booktitle_closed_string) ~isempty(idx_booktitle_open_string) ~isempty(idx_keywordlist_string) ~isempty(idx_abstract_string) ~isempty(idx_authorlist_string) ~isempty(idx_authors_string) ~isempty(idx_articletitle_string) ~isempty(idx_collectiontitle_string) ~isempty(idx_end_of_article_string) ~isempty(idx_end_of_book_article_string)];
            
            idx_string_found = find(testEverything,1);
            
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

                        [processedKeywords, Keywords, line] = getKeywords(xml,line,adj,adv,verb,noun,plnoun);

                    case 6 %% Abstract
                        
                        if isempty(idx_abstracttext_string)
                            
                            %disp('<Abstract');

                            nAbstract = nAbstract + 1;

                            [Abstract, line] = getAbstract(xml,line,adj,adv,verb,noun,plnoun);
                            
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

                        [ArticleTitle, line] = getArticleTitle(xml,line,'<ArticleTitle','</ArticleTitle',adj,adv,verb,noun,plnoun);
                        
                    case 10 %% CollectionTitle
                        
                        %disp('<CollectionTitle');
                        
                        nArticleTitle = nArticleTitle + 1;

                        [ArticleTitle, line] = getArticleTitle(xml,line,'<CollectionTitle','<CollectionTitle',adj,adv,verb,noun,plnoun);

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
            
            foundSection = false;
            
        end
        
    end
    
    %if nArticle == 2; break; end
    
end

statFileName = strcat(xmlfile(1:end-4),'-PubMed-v9.mat');

save(statFileName,'PubMed','nJournalTitle','nKeywords','nAbstract','nAuthors','nArticleTitle','-v7.3');

fclose('all');

% toc

stop = GetSecs;

disp(strcat('Time:',num2str(stop-start)));

end

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
    
end

function [newkeywords, keywords, fileline] = getKeywords(pubmedfile,fileline,adj,adv,verb,noun,plnoun)

    FID = pubmedfile;
    POSITION = ftell(FID);
    
    idx_start_of_keyword = [];
    idx_end_of_keyword = [];
    idx_end_of_section_keywordlist = [];
    
    nKeys = 1;
    
    keywords = char.empty;
    
    foundSection = false;
    while isempty(idx_end_of_section_keywordlist)
        
        idx_start_of_keyword = [];
        idx_end_of_keyword = [];
        idx_end_of_section_keywordlist = [];
    
        end_of_section_keywordlist = '</KeywordList>';
        if ~foundSection; idx_end_of_section_keywordlist = strfind(fileline,end_of_section_keywordlist); end
        if ~isempty(idx_end_of_section_keywordlist); foundSection = true; end
        
        if ~foundSection
            
            start_of_keyword = '<Keyword ';
            idx_start_of_keyword = strfind(fileline,start_of_keyword);

            end_of_keyword = '</Keyword>';
            idx_end_of_keyword = strfind(fileline,end_of_keyword);
            
            if ~isempty(idx_start_of_keyword) || ~isempty(idx_end_of_keyword); foundSection = true; end
        
        end
        
        if ~isempty(idx_start_of_keyword)
            
            end_key = '>';
            idx_end_key = strfind(fileline,end_key);
        
            keywords{nKeys} = fileline(idx_end_key(1)+1:idx_end_of_keyword-1);
        
            nKeys = nKeys + 1;
            
        end
        
        fileline = fgetl(pubmedfile);
        foundSection = false;
        
    end

    allkeywords = char.empty;
     
    if ~isempty(keywords)
        
        nKey = length(keywords);
        for iKey=1:nKey

            allkeywords = [allkeywords, ' ', keywords{iKey}];

        end
    
    end
    
    newkeywords = processParagraph(allkeywords,adj,adv,verb,noun,plnoun);
    
    if size(newkeywords,1) ~= 1; newkeywords = newkeywords'; end
    
    newkeywords = lower(newkeywords);
    
%     ORIGIN = 'bof';
%     OFFSET = POSITION + 1;
%     FID = pubmedfile;
%     STATUS = fseek(FID, OFFSET, ORIGIN);
     
end

function [abstract, fileline] = getAbstract(pubmedfile,fileline,adj,adv,verb,noun,plnoun)
        
    FID = pubmedfile;
    POSITION = ftell(FID);
    
    allAbstract = [];
    
    idx_start_of_section_abstract = [];
    idx_end_of_section_abstract = [];
    idx_start_of_sub_section_abstracttext = [];
    idx_end_of_sub_section_abstracttext = [];
    
    abstract = char.empty;
    
    foundStartSection = false;
    foundEndSection = false;
    foundStartSubSection = false;
    foundEndSubSection = false;
    openSubSection = false;
    while isempty(idx_end_of_section_abstract)
        
        idx_start_of_section_abstract = [];
        idx_end_of_section_abstract = [];
        idx_start_of_sub_section_abstracttext = [];
        idx_end_of_sub_section_abstracttext = [];
    
        start_of_section_abstract = '<Abstract';
        if ~foundStartSection; idx_start_of_section_abstract = strfind(fileline,start_of_section_abstract); end
        if ~isempty(idx_start_of_section_abstract); foundStartSection = true; end
    
        end_of_section_abstract = '</Abstract>';
        if foundStartSection; idx_end_of_section_abstract = strfind(fileline,end_of_section_abstract); end
        if ~isempty(idx_end_of_section_abstract); foundEndSection = true; end
        
        start_of_sub_section_abstracttext = '<AbstractText';
        if ~foundStartSubSection; idx_start_of_sub_section_abstracttext = strfind(fileline,start_of_sub_section_abstracttext); end
        if ~isempty(idx_start_of_sub_section_abstracttext); foundStartSubSection = true; openSubSection = true; end
        
        end_of_sub_section_abstracttext = '</AbstractText>';
        if openSubSection; idx_end_of_sub_section_abstracttext = strfind(fileline,end_of_sub_section_abstracttext); end
        if ~isempty(idx_end_of_sub_section_abstracttext); foundEndSubSection = true; end
    
        if foundStartSection
            
            start_ = idx_start_of_section_abstract+length(start_of_section_abstract)+1;
            
        end
        
        if foundStartSubSection
            
            end_string = '>';
            idx_end = strfind(fileline,end_string);
            start_ = idx_end + 1;
            
            foundStartSubSection = false;
            
        end
            
        if ~foundStartSection && ~foundStartSubSection
            
            start_ = 1;
            
        end
        
        if foundEndSection
            
            end_ = idx_end_of_section_abstract - 1; 
            
        end
        
        if foundEndSubSection
         
            end_ = idx_end_of_sub_section_abstracttext - 1;
            
        end
        
        if ~foundEndSection && ~foundEndSubSection
            
            end_ = length(fileline);
            
        end
            
        allAbstract = strcat(allAbstract,fileline(start_:end_));
        
        fileline = fgetl(pubmedfile);
        foundSection = false;

    end
    
    abstract = processParagraph(allAbstract,adj,adv,verb,noun,plnoun);

%     ORIGIN = 'bof';
%     OFFSET = POSITION + 1;
%     FID = pubmedfile;
%     STATUS = fseek(FID, OFFSET, ORIGIN);
    
end

function [authors, fileline] = getAuthors(pubmedfile,fileline,end_section) 
        
    FID = pubmedfile;
    POSITION = ftell(FID);
    
    idx_end_of_section_authorlist = [];
    idx_start_of_author_one = [];
    idx_start_of_author_two = [];
    idx_start_of_lastname = [];
    idx_end_of_lastname = [];
    idx_start_of_forename = [];
    idx_end_of_forename = [];
    
    iAuthors = 0;
    
    authors = cell.empty;
    
    foundSection = false;
    while isempty(idx_end_of_section_authorlist)
        
        idx_end_of_section_authorlist = [];
        idx_start_of_author_one = [];
        idx_start_of_author_two = [];
        idx_start_of_lastname = [];
        idx_end_of_lastname = [];
        idx_start_of_forename = [];
        idx_end_of_forename = [];
        
        idx_end_of_section_authorlist = strfind(fileline,end_section);
        if ~isempty(idx_end_of_section_authorlist); foundSection = true; end
        
        if ~foundSection; 
            
            start_of_author_one = '<Author ';
            start_of_author_two = '<Author>';
            
            idx_start_of_author_one = strfind(fileline,start_of_author_one);
        	idx_start_of_author_two = strfind(fileline,start_of_author_two);
            
            if ~isempty(idx_start_of_author_one) || ~isempty(idx_start_of_author_two); foundSection = true; end
            
        end
        
        if ~foundSection
            
            start_of_lastname = '<LastName>';
            idx_start_of_lastname = strfind(fileline,start_of_lastname);

            end_of_lastname = '</LastName>';
            idx_end_of_lastname = strfind(fileline,end_of_lastname);
            
            if ~isempty(idx_start_of_lastname) || ~isempty(idx_end_of_lastname); foundSection = true; end
        
        end
        
        if ~foundSection
            
            start_of_forename = '<ForeName>';
            idx_start_of_forename = strfind(fileline,start_of_forename);

            end_of_forename = '</ForeName>';
            idx_end_of_forename = strfind(fileline,end_of_forename);
            
            if ~isempty(idx_start_of_forename) || ~isempty(idx_end_of_forename); foundSection = true; end
        
        end
        
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
        foundSection = false;
        
    end
    
%     ORIGIN = 'bof';
%     OFFSET = POSITION + 1;
%     FID = pubmedfile;
%     STATUS = fseek(FID, OFFSET, ORIGIN);
    
end

function [articletitlekeywords, fileline] = getArticleTitle(pubmedfile,fileline,start_section,end_section,adj,adv,verb,noun,plnoun)
        
    FID = pubmedfile;
    POSITION = ftell(FID);
    
    articletitle = [];
    
    idx_start_of_section_title = [];
    idx_end_of_section_title = [];
    
    articletitlekeywords = char.empty;

    while isempty(idx_end_of_section_title)
        
        idx_start_of_section_title = [];
        idx_end_of_section_title = [];
    
        idx_start_of_section_title = strfind(fileline,start_section);
    
        idx_end_of_section_title = strfind(fileline,end_section);
    
        if ~isempty(idx_start_of_section_title)
           
            end_string = '>';
           
            idx_end = strfind(fileline,end_string);
            
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
    
    articletitlekeywords = processParagraph(articletitle,adj,adv,verb,noun,plnoun);

    
%     ORIGIN = 'bof';
%     OFFSET = POSITION + 1;
%     FID = pubmedfile;
%     STATUS = fseek(FID, OFFSET, ORIGIN);
        
end

function allkeywords = processParagraph(input_string,adj,adv,verb,noun,plnoun)
    
%     start = GetSecs;
    
    input_string = lower(input_string);
%tic
    input_string = removeUndesiredStrings(input_string);
%toc    
    input_string = strtrim(input_string);
    
    idx_spaces = isspace(input_string);
   
    allkeywords = char.empty;
%tic   
    if ~isempty(idx_spaces)
        
        input_string = strsplit(input_string);
        
%         nFirstKeys = length(input_string);
%         
%         if ~isempty(input_string)
%             
%             for i=1:nFirstKeys
% 
%                 tmp = char.empty;
%                 tmp = removeUndesiredStrings(input_string{i});
%                 tmp = strtrim(tmp);
%                 
%                 idx_spaces_tmp = isspace(tmp);
%                 
%                 if ~isempty(idx_spaces_tmp)
%                     
%                     tmp = strsplit(tmp);
%                     
%                 end
% 
%                 allkeywords = [allkeywords, tmp];
% 
%             end
%         
%         else
%             
%             allkeywords = input_string;
%             
%         end

        allkeywords = input_string;
        
    else
        
        allkeywords = input_string;
        
    end
%toc    
    nWordsInArticleTitle = length(allkeywords);
    
    idx_words_to_remove = [];
%tic    
    if isempty(idx_spaces)
        
        if ~isempty(allkeywords)
            
            idx_adj = [];
            idx_adv = [];
            idx_verb = [];
            check_if_is_plural = [];
        
            idx_adj = find(strcmpi(allkeywords,adj),1);
            idx_adv = find(strcmpi(allkeywords,adv),1);
            idx_verb = find(strcmpi(allkeywords,verb),1);
            
            if ~isempty(idx_adj) || ~isempty(idx_adv) || ~isempty(idx_verb)
            
                allkeywords = [];
            
            else
                
%                 check_if_is_plural = find(strcmpi(allkeywords,plnoun),1);
%             
%                 if ~isempty(check_if_is_plural)
%                 
%                     allkeywords = noun{check_if_is_plural};
%                 
%                 end

                 if length(allkeywors)>3
                     
                     [singleword, thereissingle, itisaplural] = checkPlural(allkeywords);

                     if itisaplural && thereissingle

                          allkeywords{iWords} = singleword;

                     elseif ~itisaplural || (itisaplural && ~thereissingle)

                          % keep the same word

                     end

                 end
                 
            end

        else
            
            allkeywords = [];
            
        end
        
    else
      
        for iWords=1:nWordsInArticleTitle
            
            idx_adj = [];
            idx_adv = [];
            idx_verb = [];
            check_if_is_plural = [];
            
%              disp('adj-adv-verb');
%              tic
            idx_adj = find(strcmpi(allkeywords{iWords},adj),1);
            idx_adv = find(strcmpi(allkeywords{iWords},adv),1);
            idx_verb = find(strcmpi(allkeywords{iWords},verb),1);
%              toc
            
            if ~isempty(idx_adj) || ~isempty(idx_adv) || ~isempty(idx_verb)

                idx_words_to_remove = [idx_words_to_remove, iWords];

            else
                
%                 disp('plural');
%                 tic
%                 
%                 check_if_is_plural = find(strcmpi(allkeywords{iWords},plnoun),1);
%                 
%                 toc
%            
%                 if ~isempty(check_if_is_plural)
%                 
%                     allkeywords{iWords} = noun{check_if_is_plural};
%                 
%                 end

                  if length(allkeywords{iWords}) > 3
     
                      [singleword, thereissingle, itisaplural] = checkPlural(allkeywords{iWords});

                      if itisaplural && thereissingle

                          allkeywords{iWords} = singleword;

                      elseif ~itisaplural || (itisaplural && ~thereissingle)

                          % keep the same word

                      end
                  
                  end
                
            end
   
        end
 
        if ~isempty(idx_words_to_remove)

            allkeywords(idx_words_to_remove) = [];

        end
        
    end
%toc      
    if size(allkeywords,1) ~= 1; allkeywords = allkeywords'; end
    
%     stop = GetSecs;
%     
%     timeDuration = stop - start;
%     
%     disp(strcat('Paragraph:',num2str(timeDuration)));


end

function output_string = removeUndesiredStrings(input_string)

    output_string = char.empty;
    
    input_string = strrep(input_string,'.',' ');
    input_string = strrep(input_string,',',' ');
    input_string = strrep(input_string,';',' ');
    input_string = strrep(input_string,'?',' ');
    input_string = strrep(input_string,'!',' ');
    input_string = strrep(input_string,' - ',' ');
    input_string = strrep(input_string,':',' ');
    
    input_string = strrep(input_string,')',' ');
    input_string = strrep(input_string,']',' ');
    input_string = strrep(input_string,'}',' ');
    input_string = strrep(input_string,'>',' ');
    input_string = strrep(input_string,'?',' ');
    
    input_string = strrep(input_string,'(',' ');
    input_string = strrep(input_string,'[',' ');
    input_string = strrep(input_string,'{',' ');
    input_string = strrep(input_string,'<',' ');
    input_string = strrep(input_string,'?',' ');
    
    input_string = strrep(input_string,'#',' ');
    input_string = strrep(input_string,'@',' ');
    input_string = strrep(input_string,'$',' ');
    input_string = strrep(input_string,'&amp',' ');
    input_string = strrep(input_string,'&quot',' ');
    input_string = strrep(input_string,'&',' ');
    input_string = strrep(input_string,'%',' ');
    
    input_string = strrep(input_string,'?',' ');
    input_string = strrep(input_string,'?',' ');
    
    input_string = strrep(input_string,'*',' ');
    input_string = strrep(input_string,'+',' ');
    input_string = strrep(input_string,'-',' ');
    input_string = strrep(input_string,'=',' ');
    input_string = strrep(input_string,'?',' ');
    input_string = strrep(input_string,'_',' ');
    input_string = strrep(input_string,'|',' ');
    input_string = strrep(input_string,'\',' ');
    input_string = strrep(input_string,'/',' ');
    input_string = strrep(input_string,'±',' ');
    input_string = strrep(input_string,'?',' ');
    input_string = strrep(input_string,'?',' ');
    
    input_string = strrep(input_string,'^',' ');
    input_string = strrep(input_string,'~',' ');
    input_string = strrep(input_string,'`',' ');
    input_string = strrep(input_string,'''',' ');
    input_string = strrep(input_string,'"',' ');
    
    input_string = strrep(input_string,'copyrightinformation',' ');
    input_string = strrep(input_string,'copyright',' ');
    input_string = strrep(input_string,'?',' ');
    
    input_string = strrep(input_string,' a ',' ');
    input_string = strrep(input_string,' an ',' ');
    input_string = strrep(input_string,' the ',' ');
    
    input_string = strrep(input_string,' than ',' ');
    input_string = strrep(input_string,' then ',' ');
    
    input_string = strrep(input_string,' this ',' ');
    input_string = strrep(input_string,' these ',' ');
    input_string = strrep(input_string,' that ',' ');
    input_string = strrep(input_string,' those ',' ');
    
    input_string = strrep(input_string,' I ',' ');
    input_string = strrep(input_string,' you ',' ');
    input_string = strrep(input_string,' he ',' ');
    input_string = strrep(input_string,' she ',' ');
    input_string = strrep(input_string,' it ',' ');
    input_string = strrep(input_string,' we ',' ');
    input_string = strrep(input_string,' they ',' ');
    
    input_string = strrep(input_string,' me ',' ');
    input_string = strrep(input_string,' my ',' ');
    input_string = strrep(input_string,' mine ',' ');
    input_string = strrep(input_string,' yours ',' ');
    input_string = strrep(input_string,' your ',' ');
    input_string = strrep(input_string,' him ',' ');
    input_string = strrep(input_string,' his ',' ');
    input_string = strrep(input_string,' her ',' ');
    input_string = strrep(input_string,' its ',' ');
    input_string = strrep(input_string,' our ',' ');
    input_string = strrep(input_string,' ours ',' ');
    input_string = strrep(input_string,' their ',' ');
    input_string = strrep(input_string,' theirs ',' ');
    input_string = strrep(input_string,' them ',' ');
    
    input_string = strrep(input_string,' and ',' ');
    input_string = strrep(input_string,' or ',' ');
    
    input_string = strrep(input_string,' into ',' ');
    input_string = strrep(input_string,' in ',' ');
    input_string = strrep(input_string,' inside ',' ');
    input_string = strrep(input_string,' out ',' ');
    input_string = strrep(input_string,' outside ',' ');
    input_string = strrep(input_string,' from ',' ');
    input_string = strrep(input_string,' to ',' ');
    input_string = strrep(input_string,' for ',' ');
    input_string = strrep(input_string,' at ',' ');
    input_string = strrep(input_string,' on ',' ');
    input_string = strrep(input_string,' above ',' ');
    input_string = strrep(input_string,' below ',' ');
    input_string = strrep(input_string,' up ',' ');
    input_string = strrep(input_string,' down ',' ');
    input_string = strrep(input_string,' upside ',' ');
    input_string = strrep(input_string,' downside ',' ');
    input_string = strrep(input_string,' underneath ',' ');
    input_string = strrep(input_string,' left ',' ');
    input_string = strrep(input_string,' right ',' ');
    input_string = strrep(input_string,' after ',' ');
    input_string = strrep(input_string,' before ',' ');
    input_string = strrep(input_string,' by ',' ');
    input_string = strrep(input_string,' of ',' ');
    
    input_string = strrep(input_string,' since ',' ');
    input_string = strrep(input_string,' until ',' ');
    input_string = strrep(input_string,' about ',' ');
    
    input_string = strrep(input_string,' what ',' ');
    input_string = strrep(input_string,' who ',' ');
    input_string = strrep(input_string,' where ',' ');
    input_string = strrep(input_string,' which ',' ');
    input_string = strrep(input_string,' whether ',' ');
    input_string = strrep(input_string,' whom ',' ');
    input_string = strrep(input_string,' whose ',' ');
    input_string = strrep(input_string,' while ',' ');
    input_string = strrep(input_string,' why ',' ');
    input_string = strrep(input_string,' when ',' ');
    input_string = strrep(input_string,' with ',' ');
    input_string = strrep(input_string,' without ',' ');
    
    input_string = strrep(input_string,' because ',' ');
    
    input_string = strrep(input_string,'0s ',' ');
    
    vowels = {' a ', ' b ', ' c ', ' d ', ' e ', ' f ', ' g ', ' h ', ' i ', ' j ', ' k ', ' l ', ' m ', ' n ', ' o ', ' p ', ' q ', ' r ', ' s ', ' t ', ' u ', ' v ', ' x ', ' y ', ' z '};
    
    numbers = {'0','1','2','3','4','5','6','7','8','9'};
    
    nVowels = length(vowels);
    
    for iVowel=1:nVowels
        
        input_string = strrep(input_string,vowels{iVowel},' ');
        
    end
    
    nNumbers = length(numbers);
    
    for iNumber=1:nNumbers
        
        input_string = strrep(input_string,numbers{iNumber},' ');
        
    end
    
    output_string = input_string;
    
    clear input_string
        
end



