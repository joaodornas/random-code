
function countString(datafile)

xml = fopen(datafile, 'r', 'l', 'UTF-8');

l = 0;
nTitle = 0;
nBookTitle = 0;
nKeywords = 0;
nAbstract = 0;
nAuthorList = 0;
nAuthors = 0;
nArticleTitle = 0;
nCollectionTitle = 0;
nPubmedArticle = 0;
nPubmedBookArticle = 0;
nBookArticleWOArticleTitle = 0;

foundBookArticle = false;
foundArticleTitle = false;

while ~feof(xml)
    
    line = fgetl(xml);
    l = l + 1;
    
    lineToBeFound = '<Title>'; %% Journal Title
    idx_title_closed_string = strfind(line,lineToBeFound);
    
    lineToBeFound = '<Title '; %% Journal Title
    idx_title_open_string = strfind(line,lineToBeFound);
    
    lineToBeFound = '<BookTitle>'; %% Journal Title
    idx_booktitle_closed_string = strfind(line,lineToBeFound);
    
    lineToBeFound = '<BookTitle '; %% Journal Title
    idx_booktitle_open_string = strfind(line,lineToBeFound);

    lineToBeFound = '<KeywordList';
    idx_keywordlist_string = strfind(line,lineToBeFound);

    lineToBeFound = '<Abstract';
    idx_abstract_closed_string = strfind(line,lineToBeFound);
    
%     lineToBeFound = '<Abstract ';
%     idx_abstract_open_string = strfind(line,lineToBeFound);

    lineToBeFound = '<AuthorList';
    idx_authorlist_string = strfind(line,lineToBeFound);
    
    lineToBeFound = '<Authors';
    idx_authors_string = strfind(line,lineToBeFound);

    lineToBeFound = '<ArticleTitle';
    idx_articletitle_string = strfind(line,lineToBeFound);
    
    lineToBeFound = '<CollectionTitle';
    idx_collectiontitle_string = strfind(line,lineToBeFound);
    
    lineToBeFound = '<PubmedArticle>';
    idx_pubmedarticle_string = strfind(line,lineToBeFound);

    lineToBeFound = '<PubmedBookArticle>';
    idx_pubmedbookarticle_string = strfind(line,lineToBeFound);
    
    lineToBeFound = '<AbstractText';
    idx_abstracttext_string = strfind(line,lineToBeFound);
    
    lineToBeFound = '</PubmedBookArticle>';
    idx_end_pubmedbookarticle_string = strfind(line,lineToBeFound);
    
    if ~isempty(idx_pubmedbookarticle_string); foundBookArticle = true; end
    
    if ~isempty(idx_articletitle_string) && foundBookArticle; foundArticleTitle = true; end
    
    if ~isempty(idx_end_pubmedbookarticle_string) && foundBookArticle && ~foundArticleTitle
        
        foundBookArticle = false;
        
        nBookArticleWOArticleTitle = nBookArticleWOArticleTitle + 1;
        
    end
       
    if ~isempty(idx_abstracttext_string); idx_abstract_closed_string = []; end

    testEverything = [~isempty(idx_title_closed_string) ~isempty(idx_title_open_string) ~isempty(idx_booktitle_closed_string) ~isempty(idx_booktitle_open_string) ~isempty(idx_keywordlist_string) ~isempty(idx_abstract_closed_string) ~isempty(idx_authorlist_string) ~isempty(idx_authors_string) ~isempty(idx_articletitle_string) ~isempty(idx_collectiontitle_string) ~isempty(idx_pubmedarticle_string) ~isempty(idx_pubmedbookarticle_string)];

    idx_string_found = find(testEverything);

    if ~isempty(idx_string_found)

        switch idx_string_found

            case 1 %% Title (Journal)

                nTitle = nTitle + 1;

            case 2 %% Title (Journal)

                nTitle = nTitle + 1;

            case 3 %% BookTitle (Journal)

                nBookTitle = nBookTitle + 1;

            case 4 %% BookTitle (Journal)

                nBookTitle = nBookTitle + 1;

            case 5 %% KeywordList

                nKeywords = nKeywords + 1;

            case 6 %% Abstract
                
                if isempty(idx_abstracttext_string)

                    nAbstract = nAbstract + 1;
                    
                end
            
%             case 7 %% Abstract
% 
%                 nAbstract = nAbstract + 1;

            case 7 %% AuthorList

                nAuthorList = nAuthorList + 1;
                
            case 8 %% Authors

                nAuthors = nAuthors + 1;

            case 9 %% ArticleTitle

                nArticleTitle = nArticleTitle + 1;       
                
            case 10 %% CollectionTitle

                nCollectionTitle = nCollectionTitle + 1;
                
            case 11 %% PubmedArticle
                
                nPubmedArticle = nPubmedArticle + 1;
                
            case 12 %% PubmedBookArticle
                
                nPubmedBookArticle = nPubmedBookArticle + 1;
                
        end
        
    end
    
end

fclose('all');

disp(strcat('Number of lines:',int2str(l)));
disp(strcat('Number of nTitle:',int2str(nTitle)));
disp(strcat('Number of nBookTitle:',int2str(nBookTitle)));
disp(strcat('Number of nKeywords:',int2str(nKeywords)));
disp(strcat('Number of nAbstract:',int2str(nAbstract)));
disp(strcat('Number of nAuthorList:',int2str(nAuthorList)));
disp(strcat('Number of nAuthors:',int2str(nAuthors)));
disp(strcat('Number of nArticleTitle:',int2str(nArticleTitle)));
disp(strcat('Number of nCollectionTitle:',int2str(nCollectionTitle)));
disp(strcat('Number of nPubmedArticle:',int2str(nPubmedArticle)));
disp(strcat('Number of nPubmedBookArticle:',int2str(nPubmedBookArticle)));
disp(strcat('Number of nBookArticleWOArticleTitle:',int2str(nBookArticleWOArticleTitle)));

end