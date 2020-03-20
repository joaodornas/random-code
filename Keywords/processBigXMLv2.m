function processBigXMLv2(xmlfile)

xml = fopen(xmlfile, 'r', 'l', 'UTF-8');

fileToWrite = fopen(strcat(xmlfile(1:end-4),'-keywords.xml'), 'w', 'l', 'UTF-8');

expectedLine{1} = '<?xml version="1.0"?>';
    
expectedLine{2} = '<PubmedArticleSet>';
    
expectedLine{3} = '<PubmedArticle>';
    
expectedLine{4} = '<KeywordList';
     
expectedLine{5} = '<Keyword';

expectedLine{6} = '</Keyword>';
      
expectedLine{7} = '</KeywordList>';
      
expectedLine{8} = '</PubmedArticle>';
      
expectedLine{9} = '</PubmedArticleSet>';

l = 1;

fwrite(fileToWrite, expectedLine{1}, 'char', 0, 's');
fprintf(fileToWrite, '\n');

fwrite(fileToWrite, expectedLine{2}, 'char', 0, 's');
fprintf(fileToWrite, '\n');

ksearch = 0;

nArticle = 0;
nKeywordList = 0;

while ~feof(xml)
    
    if l == 1; line = fgetl(xml); end
    
    if mod(l,100000) == 0; disp(strcat('start loop line:',int2str(l/100000))); end
    
    idx_keywordlist_string = [];
    
    %disp('Searching for keyword list');
    
    %tic
    
    ksearch = ksearch + 1;
    
    while isempty(idx_keywordlist_string)
    
        iExpectedLine = 3;
        idx_article_string = strfind(line,expectedLine{iExpectedLine});
        
        if ~isempty(idx_article_string); nArticle = nArticle + 1; end
        
        iExpectedLine = 4;
        idx_keywordlist_string = strfind(line,expectedLine{iExpectedLine});
        
        if ~isempty(idx_keywordlist_string); 
            
            nKeywordList = nKeywordList + 1; 
            keywordListLine = line;
            
        end
        
        line = fgetl(xml);
        
        l = l + 1;
        
        if mod(l,100000) == 0; disp(strcat('search keyword list line:',int2str(l/100000))); end
        
        %disp(strcat('ksearch:',int2str(ksearch)));
        
        %position = ftell(xml);
        
        %disp(strcat('position:',int2str(position)));
        
        if feof(xml); break; end
    
    end
    
    %toc
    
    iExpectedLine = 3;
    fwrite(fileToWrite, expectedLine{iExpectedLine}, 'char', 0, 's');
    fprintf(fileToWrite, '\n');
    
    iExpectedLine = 4;
    fwrite(fileToWrite, keywordListLine, 'char', 0, 's');
    fprintf(fileToWrite, '\n');
    
    isTheEndOfKeywordList = false;
    
    %disp('Reading keyword list');
    
    %tic
    
    while ~isTheEndOfKeywordList

        iExpectedLine = 5;
        idx_keyword_string = strfind(line,expectedLine{iExpectedLine});
        
        iExpectedLine = 7;
        idx_keywordlist_string = strfind(line,expectedLine{iExpectedLine});

        if ~isempty(idx_keyword_string)
            
            keyword = line(idx_keyword_string:end);

            fwrite(fileToWrite, keyword, 'char', 0, 's');
            fprintf(fileToWrite, '\n');
            
        end
        
        if ~isempty(idx_keywordlist_string)
            
            iExpectedLine = 7;
            fwrite(fileToWrite, expectedLine{iExpectedLine}, 'char', 0, 's');
            fprintf(fileToWrite, '\n');
            
            iExpectedLine = 8;
            fwrite(fileToWrite, expectedLine{iExpectedLine}, 'char', 0, 's');
            fprintf(fileToWrite, '\n');
            
            isTheEndOfKeywordList = true;
            
        end
        
        line = fgetl(xml);
        
        l = l + 1;
        
        if mod(l,100000) == 0; disp(strcat('reading keyword list line:',int2str(l/100000))); end
        
        if feof(xml); break; end
        
    end
    
    %toc
    
end

fwrite(fileToWrite, expectedLine{9}, 'char', 0, 's');
fprintf(fileToWrite, '\n');

statFileName = strcat(xmlfile(1:end-4),'-nArticle-nKeywordList.mat');

save(statFileName,'nArticle','nKeywordList');

fclose(fileToWrite);

fclose('all');

end

