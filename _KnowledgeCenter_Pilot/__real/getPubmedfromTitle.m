
baseSearchURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?';
baseFetchURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';

time_to_pause = 60;

%%% MENDELEY
% tFile_int = [1000 2000 3000 4000 5000 6000 7000];
% tFile_int = [3000 4000 5000 6000 7000];
tFile_int = [4000 5000 6000 7000];

for iMen=tFile_int
    
   load(strcat('Mendeley_',int2str(iMen),'.mat')); 
   
   nTitles = length(title);
   
   if iMen == 4000; start_title = 600; else start_title = 1; end
   
   for iTitle=start_title:nTitles
       
       this_title = title{iTitle};
       
       disp(strcat('Mendeley_',int2str(iMen),'-',int2str(iTitle)));
       
       if this_title ~= 0 & ~isempty(file{iTitle})
          
           this_title = strrep(this_title,'"','');
           this_title = strrep(this_title,' ','+');
           this_title = strrep(this_title,',','');
           this_title = lower(this_title);
           
           dbsource = 'db=pubmed';
           term = strcat('&term=',this_title);
           field = '&field=title';
           usehistory = '&usehistory=y';
           retstart = '&retstart=0';
           retmax = '&retmax=200';
           rettype = '&rettype=uilist';
           datetype = '&datetype=pdat';
           
           mindate = strcat('&mindate=','1900');
           maxdate = strcat('&maxdate=','2016');
           
           noFinished = true;help 
           
           while noFinished
               
               try 

                       searchURL = [baseSearchURL,dbsource,term,field,usehistory,retstart,retmax,rettype,datetype,mindate,maxdate];
                       [whole_html, status_search] = urlread(searchURL);

                       ncbi = regexp(whole_html,'<QueryKey>(?<QueryKey>\w+)</QueryKey>\s*<WebEnv>(?<WebEnv>\S+)</WebEnv>','names');

                       webenv = strcat('&webenv=',ncbi.WebEnv); 
                       querykey = strcat('&query_key=',ncbi.QueryKey); 

                       rettype = '&rettype=';
                       retmode = '&retmode=xml';

                       fetchURL = [baseFetchURL,dbsource,webenv,querykey,term,field,usehistory,retstart,retmax,rettype,retmode,datetype,mindate,maxdate];

                       [whole_html, status_fetch] = urlread(fetchURL);
                       
                       noFinished = false;

               catch
                   
                   disp('Error on Internet connection. Retrying in 60 seconds.');
                   
                   pause(time_to_pause);

               end
           
           end
           
           this_file = file(iTitle);
           idx_slash = strfind(this_file{:},'/');
           tmp = this_file{:};
           this_file = tmp(idx_slash(end)+1:end);
           
           this_file = strrep(this_file,'"','');
           this_file = strrep(this_file,'''','');
           
           if status_search & status_fetch
               
                fid = fopen(strcat(this_file(1:end-4),'.html'),'w');
                fprintf(fid,'%s\n',whole_html);
                fclose(fid);
                
           else
               
              fid = fopen(strcat(this_file(1:end-4),'.html'),'w');
              fprintf(fid,'%s\n',strcat('ErrorStatus:',int2str(status_search),':',int2str(status_fetch)));
              fclose(fid);
              
              disp(strcat('ErrorStatus:',int2str(status_search),':',int2str(status_fetch)));
               
           end
           
           pause(time_to_pause);
           
       end
       
   end
    
end
