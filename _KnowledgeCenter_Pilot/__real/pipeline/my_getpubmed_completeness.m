function my_getpubmed_completeness

% Intervals between downloads
minutes_to_pause = 5;
seconds_to_pause = 60*minutes_to_pause;
time_between_periods = 24*60*60;

% Create base URL for PubMed db site 
baseSearchURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?';
baseFetchURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';

% Generate Dates
% % start_year = 1900;
% % start_month = 01;
% % start_day = 01;

start_year = 1924;
start_month = 07;
start_day = 12;

end_year = 2017;
end_month = 12;
end_day = 31;

first_day = datetime(start_year,start_month,start_day);
last_day = datetime(end_year,end_month,end_day);

period = first_day:last_day;

total_days = length(period);

my_dates = zeros(total_days,3);

for iD=1:total_days
    
    this_day = period(iD);
    
    my_dates(iD,1) = this_day.Day;
    my_dates(iD,2) = this_day.Month;
    my_dates(iD,3) = this_day.Year;
    
end

nDates = length(my_dates);

% Load search strings
dbsource = 'db=pubmed';
term = '&term=';
field = '&field=';
usehistory = '&usehistory=y';
retstart = '&retstart=0';
retmax = '&retmax=200';
rettype = '&rettype=uilist';
datetype = '&datetype=pdat';

iSteps = 200;

for iDate=1:nDates

    tic;
    
        start_day = my_dates(iDate,1);
        start_month = my_dates(iDate,2);
        start_year = my_dates(iDate,3);
  
        start_day_string = strcat('0',int2str(start_day));
        end_day_string = strcat('0',int2str(start_day));

        if start_month < 10

            month_string = strcat('0',int2str(start_month));

        else

            month_string = int2str(start_month);

        end

        start_this_day_search = strcat(int2str(start_year),'/',month_string,'/',start_day_string);
        end_this_day_search = strcat(int2str(start_year),'/',month_string,'/',end_day_string);

        mindate = strcat('&mindate=',start_this_day_search);
        maxdate = strcat('&maxdate=',end_this_day_search);

        % Create search URL
        searchURL = [baseSearchURL,dbsource,term,field,usehistory,retstart,retmax,rettype,datetype,mindate,maxdate];

        done = 0;
        while ~done
            
            try
                [medlineText STATUS] = urlread(searchURL);
                
                if STATUS
                    
                    done = 1;
                    
                end
                
                disp(strcat('pausing 5 minutes:',datestr(now)));
                pause(seconds_to_pause);
                
            catch ME
                getReport(ME);
                disp(ME);
                disp(strcat('error,pausing:',datestr(now)));
                pause(seconds_to_pause);
            end
            
        end

        counts = regexp(medlineText, '<Count>(\w*)</Count>', 'tokens');
        
        if ~isempty(counts)
            
           nTotalPapers = str2num(cell2mat(counts{1}));
            
           if nTotalPapers > iSteps
               
               nSteps = (nTotalPapers - mod(nTotalPapers,iSteps))/iSteps;
               
           else
               
               nSteps = 0;
               
           end
           
           allSteps = (0:nSteps).*iSteps;
           
           if allSteps(end) == nTotalPapers; allSteps(end) = []; end

            ncbi = regexp(medlineText,'<QueryKey>(?<QueryKey>\w+)</QueryKey>\s*<WebEnv>(?<WebEnv>\S+)</WebEnv>','names');

            webenv = strcat('&webenv=',ncbi.WebEnv); 
            querykey = strcat('&query_key=',ncbi.QueryKey); 

            rettype = '&rettype=';
            retmode = '&retmode=xml';

            for iRetStart=allSteps

                retstart = strcat('&retstart=',int2str(iRetStart));

                fetchURL = [baseFetchURL,dbsource,webenv,querykey,term,field,usehistory,retstart,retmax,rettype,retmode,datetype,mindate,maxdate];

                
                done = 0;
                while ~done
            
                    try
                        [medlineText, STATUS] = urlread(fetchURL);
                        
                        if STATUS
                            
                            done = 1;
                            
                        end
                        
%                         disp(strcat('pausing:',datestr(now)));
%                         pause(seconds_to_pause);
                        
                    catch ME
                        getReport(ME);
                        disp(ME);
                        disp(strcat('error,pausing:',datestr(now)));
                        pause(seconds_to_pause);
                    end
            
                end
        
                disp(strcat('saving, pausing:',datestr(now)));
                pause(seconds_to_pause); %% PAUSE a little bit

                % medlineText = medlineText(31<medlineText &
                % medlineText<127); This line removes ACCENTUATION
                
                medlineText = medlineText(medlineText>31); % This line
                % removes CONTROL CHARACTERS

                fileID = fopen(strcat(int2str(start_year),'-',month_string,'-',start_day_string,'-',int2str(iRetStart),'.xml'),'w','native','utf-8');
                nbytes = fprintf(fileID,'%s',medlineText);
                fclose(fileID);

                % [tree, RootName, DOMnode] = xml_read(strcat(label,'-',int2str(query_id),'-',int2str(iYear),'-',month_string,'-',int2str(iRetStart),'.xml'));
                % save(strcat(label,'-',int2str(query_id),'-',int2str(iYear),'-',month_string,'-',int2str(iRetStart),'.mat'),'tree');

            end
            
        end
        
        t = toc;
        
        if t > time_between_periods
            
            tic;
            
            disp(strcat('pausing 1 day:',datestr(now)));
            pause(time_between_periods);
            
        end

    end

end
   

