function processBigXMLv6(xmlfile)

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

wordNet = load('noun.mat');
noun = lower(strtrim(wordNet.uniqueWords));

wordNet = load('noun-pl.mat');
plnoun = strtrim(wordNet.plNouns);

nArticle = 0;
nJournalTitle = 0;
nKeywords = 0;
nAuthors = 0;
nArticleTitle = 0;
nAbstract = 0;

ArticleTitle = char.empty;
Abstract = char.empty;
Keywords = char.empty;

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
        
        if nArticle > 1
            
            allKeywords = [ArticleTitle, Abstract, Keywords];
      
            if ~isempty(allKeywords)
            
                [allKeywords, nallKeywords] = count_unique(allKeywords);
                
            end
            
            PubMed(nArticle - 1).allKeywords = allKeywords;
            
            allKeywords = char.empty;
            ArticleTitle = char.empty;
            Abstract = char.empty;
            Keywords = char.empty;
            
        end
            
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

                        [Keywords, line] = getKeywords(xml,line);

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
            
        end
        
    end
    
end

statFileName = strcat(xmlfile(1:end-4),'-PubMed-v6.mat');

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

function [abstract, fileline] = getAbstract(pubmedfile,fileline,adj,adv,verb,noun,plnoun)
        
    FID = pubmedfile;
    POSITION = ftell(FID);
    
    allAbstract = [];
    idx_start_of_section_abstract = [];
    idx_end_of_section_abstract = [];
    
    abstract = char.empty;
    
    while isempty(idx_end_of_section_abstract)
    
        start_of_section_abstract = '<Abstract';
        idx_start_of_section_abstract = strfind(fileline,start_of_section_abstract);
    
        end_of_section_abstract = '</Abstract>';
        idx_end_of_section_abstract = strfind(fileline,end_of_section_abstract);
        
        start_of_sub_section_abstracttext = '<AbstractText';
        idx_start_of_sub_section_abstracttext = strfind(fileline,start_of_sub_section_abstracttext);
        
        end_of_sub_section_abstracttext = '</AbstractText>';
        idx_end_of_sub_section_abstracttext = strfind(fileline,end_of_sub_section_abstracttext);
    
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
        
        if ~isempty(idx_end_of_sub_section_abstracttext)
            
            end_ = idx_end_of_sub_section_abstracttext - 1;
            
        else
            
            end_ = length(fileline);
            
        end
        
        if ~isempty(idx_end_of_section_abstract); 
        
            end_ = idx_end_of_section_abstract - 1; 
    
        else
        
            end_ = length(fileline);
    
        end
            
        allAbstract = strcat(allAbstract,fileline(start_:end_));
        
        fileline = fgetl(pubmedfile);
        
    end
    
    abstract = processParagraph(allAbstract,adj,adv,verb,noun,plnoun);

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

function [articletitlekeywords, fileline] = getArticleTitle(pubmedfile,fileline,start_section,end_section,adj,adv,verb,noun,plnoun)
        
    FID = pubmedfile;
    POSITION = ftell(FID);
    
    articletitle = [];
    idx_start_of_section_title = [];
    idx_end_of_section_title = [];
    
    articletitlekeywords = char.empty;
    
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
    
    articletitlekeywords = processParagraph(articletitle,adj,adv,verb,noun,plnoun);

    
%     ORIGIN = 'bof';
%     OFFSET = POSITION + 1;
%     FID = pubmedfile;
%     STATUS = fseek(FID, OFFSET, ORIGIN);
        
return

function allkeywords = processParagraph(input_string,adj,adv,verb,noun,plnoun)

    start = GetSecs;
    
    input_string = lower(input_string);

    input_string = removeUndesiredStrings(input_string);
   
    idx_spaces = find(isspace(input_string));
    
    selectedKeys = char.empty;
    allkeywords = char.empty;
    
    if ~isempty(idx_spaces)
        
        selectedKeys{1} = input_string(1:idx_spaces(1)-1);
    
        for i=1:length(idx_spaces)-1
       
            selectedKeys{i+1} = input_string(idx_spaces(i)+1:idx_spaces(i+1)-1);
        
        end
    
        selectedKeys{length(idx_spaces)+1} = input_string(idx_spaces(end)+1:end);
    
        allkeywords = selectedKeys;

    else
        
        allkeywords{1} = input_string;
        
    end
    
    idx_plnoun = [];
    if ~isempty(allkeywords)
        
        if isempty(idx_spaces)
            
            idx_plnoun = find(strcmpi(allkeywords,plnoun));
            
            if ~isempty(idx_plnoun)
                
                if length(allkeywords) > 3
                    
                    idx_noun = find(strmatch(allkeywords(1:end-3),noun));
                
                    if ~isempty(idx_noun)
                   
                        allkeywords{1} = noun{idx_noun};
                        nallkeywords = 1;
                   
                    else
                   
                        allkeywords{1} = plnoun{idx_plnoun};
                        nallkeywords = 1;
                   
                    end 
                    
                else
                    
                    allkeywords{1} = plnoun{idx_plnoun};
                    nallkeywords = 1;
                        
                end
                
            end
            
        else
            
            for iWords=1:length(allkeywords)
                
               idx_plnoun = find(strcmpi(allkeywords{iWords},plnoun));
                
               if ~isempty(idx_plnoun)
                   
                   if length(allkeywords{iWords}) > 3
                   
                        idx_noun = find(strmatch(allkeywords{iWords}(1:end-3),noun));
                
                        if ~isempty(idx_noun)
                        
                            allkeywords{iWords} = noun{idx_noun};
                        
                        else
                        
                            allkeywords{iWords} = plnoun{idx_plnoun};
                        
                        end 
                        
                   else
                       
                       allkeywords{iWords} = plnoun{idx_plnoun};
                       
                   end
                
               end
               
            end
        
            [allkeywords, nallkeywords] = count_unique(allkeywords);
            
        end
    
    end
    
    nWordsInArticleTitle = length(allkeywords);
    
    idx_adj = [];
    idx_adv = [];
    idx_verb = [];
    
    idx_words_to_remove = [];
    
    if isempty(idx_spaces)
        
        if ~isempty(allkeywords)
        
            idx_adj = find(strmatch(allkeywords(1:end-2),adj));
            idx_adv = find(strmatch(allkeywords(1:end-2),adv));
            idx_verb = find(strmatch(allkeywords(1:end-2),verb));
            
            if ~isempty(idx_adj) || ~isempty(idx_adv) || ~isempty(idx_verb)
            
                allkeywords = [];
            
            end

        else
            
            allkeywords = [];
            
        end
        
    else
        
        for iWords=1:nWordsInArticleTitle

            idx_adj = find(strmatch(allkeywords{iWords}(1:end-2),adj));
            idx_adv = find(strmatch(allkeywords{iWords}(1:end-2),adv));
            idx_verb = find(strmatch(allkeywords{iWords}(1:end-2),verb));

            if ~isempty(idx_adj) || ~isempty(idx_adv) || ~isempty(idx_verb)

                idx_words_to_remove = [idx_words_to_remove, iWords];

            end
   
        end

        if ~isempty(idx_words_to_remove)

            allkeywords(idx_words_to_remove) = [];

        end
        
    end
    
    if size(allkeywords,1) ~= 1; allkeywords = allkeywords'; end
    
    stop = GetSecs - start;
    
    %disp(strcat('Paragraph took:',num2str(stop),'seconds'));


return

function output_string = removeUndesiredStrings(input_string)

    output_string = char.empty;
    
    input_string = strrep(input_string,'.',' ');
    input_string = strrep(input_string,',',' ');
    input_string = strrep(input_string,';',' ');
    input_string = strrep(input_string,'?',' ');
    input_string = strrep(input_string,'!',' ');
    input_string = strrep(input_string,'-',' ');
    
    input_string = strrep(input_string,')',' ');
    input_string = strrep(input_string,']',' ');
    input_string = strrep(input_string,'}',' ');
    
    input_string = strrep(input_string,'(',' ');
    input_string = strrep(input_string,'[',' ');
    input_string = strrep(input_string,'{',' ');
    input_string = strrep(input_string,'<',' ');
    input_string = strrep(input_string,'>',' ');
    
    input_string = strrep(input_string,':','');
    input_string = strrep(input_string,'-','');
    
    input_string = strrep(input_string,'#',' ');
    input_string = strrep(input_string,'@',' ');
    input_string = strrep(input_string,'$',' ');
    input_string = strrep(input_string,'&amp',' ');
    input_string = strrep(input_string,'&quot',' ');
    input_string = strrep(input_string,'&',' ');
    
    input_string = strrep(input_string,'*',' ');
    input_string = strrep(input_string,'+',' ');
    input_string = strrep(input_string,'=',' ');
    input_string = strrep(input_string,'_',' ');
    input_string = strrep(input_string,'"',' ');
    input_string = strrep(input_string,'|',' ');
    input_string = strrep(input_string,'\',' ');
    input_string = strrep(input_string,'/',' ');
    
    input_string = strrep(input_string,'^',' ');
    input_string = strrep(input_string,'~',' ');
    input_string = strrep(input_string,'`',' ');
    input_string = strrep(input_string,'''',' ');
    input_string = strrep(input_string,'?',' ');
    
    input_string = strrep(input_string,'?',' ');
    input_string = strrep(input_string,'copyrightinformation',' ');
    input_string = strrep(input_string,'copyright',' ');
    
    
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
        
return



