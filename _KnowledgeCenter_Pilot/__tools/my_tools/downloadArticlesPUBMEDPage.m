
downloadPath = 'D:\__PROJECTS\__Retina\_ArticleTitle_v4\__HTML-MISSING-PDF';
filetypes = {'pdf'};

seconds_to_pause = 60;

for iArticle=1:length(articles);
    
    fetchURL = articles{iArticle,2};
    
    if ~strcmp('NoField',fetchURL)
    
        PMID = str2num(articles{iArticle,1});

        done = 0;
        while ~done

            try
                
                website = fetchURL;
                crawling(website, filetypes, downloadPath);
                
                % [medlineText, STATUS] = urlread(fetchURL);

%                 if STATUS
% 
%                     done = 1;
% 
%                     fid = fopen(strcat(int2str(PMID),'.html'),'w');
%                     fprintf(fid,'%s\n',medlineText);
%                     fclose(fid);
% 
%                 end

                done = 1;

            catch ME

                getReport(ME);
                disp(ME);

            end

            pause(seconds_to_pause);

        end
        
    end
    
end