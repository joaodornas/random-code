
function whatDoesItHaveWO(datafile,searchedString)

xmlfile = fopen(datafile, 'r', 'l', 'UTF-8');

fileToWrite = fopen(strcat(xmlfile(1:end-4),'-without-',searchedString{1},'.xml'), 'w', 'l', 'UTF-8');

pubmedArticle = char.empty;
iLine = 0;
l = 0;

insideArticle = false;
foundSearchedString = false;

while ~feof(xmlfile)
    
    line = fgetl(xmlfile);
    l = l + 1;
    
    idx_searched_string = [];
    
    for s=1:length(searchedString)
        
        idx = strfind(line,searchedString{s});
        
        if ~isempty(idx); 
        
            idx_searched_string = idx;
        
        end
        
    end
    
    lineToBeFound = '<PubmedArticle>';
    idx_pubmedarticle_open_string = strfind(line,lineToBeFound);

    lineToBeFound = '<PubmedBookArticle>';
    idx_pubmedbookarticle_open_string = strfind(line,lineToBeFound);
    
    lineToBeFound = '</PubmedArticle>';
    idx_pubmedarticle_close_string = strfind(line,lineToBeFound);

    lineToBeFound = '</PubmedBookArticle>';
    idx_pubmedbookarticle_close_string = strfind(line,lineToBeFound);
    
    testEverything = [~isempty(idx_searched_string) ~isempty(idx_pubmedarticle_open_string) ~isempty(idx_pubmedbookarticle_open_string) ~isempty(idx_pubmedarticle_close_string) ~isempty(idx_pubmedbookarticle_close_string)];

    idx_string_found = find(testEverything);

    if ~isempty(idx_string_found)

        switch idx_string_found

            case 1 %% SearchedString
    
                 iLine = iLine + 1;
                 
                 pubmedArticle{iLine} = line;
                 
                 foundSearchedString = true;
                

            case 2 %% PubmedArticle
                
                iLine = iLine + 1;

                pubmedArticle{iLine} = line;
                
                insideArticle = true;

            case 3 %% PubmedBookArticle

                iLine = iLine + 1;

                pubmedArticle{iLine} = line;
                
                insideArticle = true;
                

            case 4 %% /PubmedArticle

                iLine = iLine + 1;
                
                pubmedArticle{iLine} = line;
                
                insideArticle = false;
                
                if foundSearchedString
                    
                    pubmedArticle = char.empty;
                    
                    foundSearchedString = false;
                    
                    iLine = 0;
                    
                else
                    
                    for i=1:iLine
                        
                        fwrite(fileToWrite,pubmedArticle{i}, 'char', 0, 's');
                        fprintf(fileToWrite, '\n');
                        
                    end
                    
                    return
                    
                end
                
            case 5 %% /PubmedBookArticle

                iLine = iLine + 1;
                
                pubmedArticle{iLine} = line;
                
                insideArticle = false;
                
                if foundSearchedString
                    
                    pubmedArticle = char.empty;
                    
                    foundSearchedString = false;
                    
                    iLine = 0;
                    
                else
                    
                    for i=1:iLine
                        
                        fwrite(fileToWrite,pubmedArticle{i}, 'char', 0, 's');
                        fprintf(fileToWrite, '\n');
                        
                    end
                    
                    return
                    
                end
                
        end
        
    elseif insideArticle
       
        iLine = iLine + 1;
        
        pubmedArticle{iLine} = line;
        
    end
    
end

fclose('all');

end