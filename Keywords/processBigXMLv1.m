function processBigXMLv1(xmlfile)


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

while ~feof(xml)
    
    if mod(l,100000) == 0; disp(strcat('line:',int2str(l/100000))); end
    
    line = fgetl(xml);
    
    foundKeywordListNow = false;
    
    for iExpectedLine=1:length(expectedLine)
        
        lengthExpectedLine = length(expectedLine{iExpectedLine});
        
        idx_start_string = strfind(line,expectedLine{iExpectedLine});
        
        if ~isempty(idx_start_string)
            
            if iExpectedLine == 4
                
                expectedLineToWrite = line;
                foundKeywordListNow = true;
                
            else
                
                expectedLineToWrite = expectedLine{iExpectedLine};
                
            end
            
            if ~(foundKeywordListNow && iExpectedLine == 5)
                
                fwrite(fileToWrite, expectedLineToWrite, 'char', 0, 's');
        
            end
            
            if iExpectedLine == 5
                
                if ~foundKeywordListNow
                    
                    idx_end_string = strfind(line,expectedLine{6});
                
                    if isempty(idx_end_string); idx_end_string = length(line) + 1; end
                
                    keyword = line(idx_start_string+lengthExpectedLine:idx_end_string-1);
                
                    fwrite(fileToWrite,keyword, 'char', 0, 's');
                    
                else
                    
                    foundKeywordListNow = false;
                    
                end
                
            else
                
                fprintf(fileToWrite, '\n');
                
            end
            
        end
           
    end
    
    l = l + 1;
    
end

l - 1

fclose(fileToWrite);

%fclose(xmlfile);

fclose('all');

end

