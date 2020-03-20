% function my_getpubmed
% 
% query_id = 1;
% 
% start_year = 1900;
% end_year = 2016;
% getpubmed('retina',query_id,start_year,end_year);
% 
% start_year = 1998;
% end_year = 2016;
% 
% getpubmed('visual',query_id,start_year,end_year);
% 
% start_year = 1900;
% end_year = 2016;
% 
% getpubmed('attention',query_id,start_year,end_year);
% getpubmed('awareness',query_id,start_year,end_year);
% getpubmed('consciousness',query_id,start_year,end_year);
% 
% end

function my_getpubmed(label,query_id,start_year,end_year,minutes_to_pause)

% Create base URL for PubMed db site 
baseSearchURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?';
baseFetchURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';

%%% VISUAL SYSTEM
% if strcmp(label,'retina'); searchterm = '(retina) OR (retina AND photoreceptor) OR (retina AND cone) OR (retina AND rod) OR (retina AND bipolar) OR (retina AND horizontal) OR (retina AND amacrine) OR (retina AND ganglion)'; end
if strcmp(label,'retina'); searchterm = '(retina) OR (retina AND photoreceptor) OR (retina AND cone) OR (retina AND rod) OR (retina AND bipolar) OR (retina AND horizontal) OR (retina AND amacrine) OR (retina AND ganglion)'; end
if strcmp(label,'retinal'); searchterm = '(retinal) OR (retinal AND photoreceptor) OR (retinal AND cone) OR (retinal AND rod) OR (retinal AND bipolar) OR (retinal AND horizontal) OR (retinal AND amacrine) OR (retinal AND ganglion)'; end

if strcmp(label,'visual'); searchterm = 'visual'; end

%%% AUDITORY SYSTEM
if strcmp(label,'cochlea?'); searchterm = '???'; end
if strcmp(label,'auditory'); searchterm = 'auditory'; end

%%% MOTOR SYSTEM
if strcmp(label,'motor'); searchterm = 'motor'; end

%%% FUNCTIONS
if strcmp(label,'memory'); searchterm = 'memory'; end
if strcmp(label,'learning'); searchterm = 'learning'; end
if strcmp(label,'decision_making'); searchterm = 'decision AND making'; end
if strcmp(label,'emotions'); searchterm = 'emotions OR emotion'; end

%%% STATES
if strcmp(label,'sleep'); searchterm = 'sleep'; end
if strcmp(label,'awake'); searchterm = 'awake'; end
if strcmp(label,'resting_state'); searchterm = 'resting AND state'; end
if strcmp(label,'attention'); searchterm = 'attention'; end

%%% MIND NEUROSCIENCE
if strcmp(label,'awareness'); searchterm = 'awareness'; end
if strcmp(label,'consciousness'); searchterm = 'consciousness'; end

%%% NEUROTRANSMITTERS
% if strcmp(label,'neurotransmitters'); searchterm = 'neurotransmitter OR endorphin OR oxytocin OR octopamine OR nicotine OR neuropharmacology OR neuropeptide OR nephrin OR muscarine OR morphine OR epinephrine OR endorphin OR histamine OR histidine OR glycerol OR glutamine OR glutamate OR epinephrine OR epinephrine OR endocannabinoid OR desipramine OR cholinesterase OR cannabinoid OR calmodulin OR caffeine OR benzodiazepine OR angiotensin OR amphetamine OR adenosine OR synapses OR synapse OR adenosylmethionine OR oxytocin OR muscarine OR histidine OR histamine OR endorphin OR endocannabinoid OR catecholamine OR norepinephrine OR catecholamine OR benzodiazepine OR barbiturates OR aspartate OR arginine OR antagonist OR agonist OR amphetamine OR adenosine OR adenine OR kainate OR noradrenaline OR acth OR acetyl OR acetylcholine OR epinephrine OR dopa OR dopamine OR dopamine OR gaba OR gabaa OR gabab OR gabaa OR methyl OR nmda OR gly OR glycine OR adrenaline'; end
if strcmp(label,'neurotransmitters'); searchterm = 'neurotransmitters OR neurotransmitter'; end

%%% SIGNAL + METHOD
if strcmp(label,'eeg'); searchterm = 'eeg OR electroencephalography'; end
if strcmp(label,'meg'); searchterm = 'meg OR magnetoencephalography'; end
if strcmp(label,'fmri'); searchterm = 'fmri OR DTI OR functional magnetic resonance imaging OR diffusion tractography imaging'; end
if strcmp(label,'optogenetic'); searchterm = 'optogenetic OR optogenetics'; end
if strcmp(label,'electrophysiology'); searchterm = 'electrophysiology OR single unit OR multi unit'; end
if strcmp(label,'microscopy'); searchterm = 'microscopy'; end
if strcmp(label,'clarity'); searchterm = 'clarity'; end

%%% FUNCTIONAL NETWORKS
if strcmp(label,'DAN'); searchterm = 'DAN OR (Dorsal AND Attention AND Network)'; end
if strcmp(label,'VAN'); searchterm = 'VAN OR (Ventral AND Attention AND Network)'; end
if strcmp(label,'VIS'); searchterm = 'VIS OR (Visual AND Network)'; end
if strcmp(label,'AUD'); searchterm = 'AUD OR (Auditory AND Network)'; end
if strcmp(label,'LAN'); searchterm = 'LAN OR (Language AND Network)'; end
if strcmp(label,'FPC'); searchterm = 'FPC OR (Frontal AND Parietal AND Control AND Network)'; end
if strcmp(label,'SMN'); searchterm = 'SMN OR (Sensory AND Motor AND Network)'; end
if strcmp(label,'DMN'); searchterm = 'DMN OR (Default AND Mode AND Network)'; end

%%% AAL
if strmatch('AAL',label); searchterm = strrep(label(5:end-2),'_',' AND '); end

disp(searchterm);

dbsource = 'db=pubmed';
term = strcat('&term=',searchterm);
field = '&field=';
usehistory = '&usehistory=y';
retstart = '&retstart=0';
retmax = '&retmax=200';
rettype = '&rettype=uilist';
datetype = '&datetype=pdat';

% start_year = 1900;
% end_year = 2016;
start_day = 1;
end_day = 31;
start_month = 1;
end_month = 12;

seconds_to_pause = 60*minutes_to_pause;

iSteps = 200;

for iYear=start_year:end_year

    for iMonth=start_month:end_month

        start_day_string = strcat('0',int2str(start_day));
        end_day_string = int2str(end_day);

        if iMonth < 10

            month_string = strcat('0',int2str(iMonth));

        else

            month_string = int2str(iMonth);

        end

        start_this_day_search = strcat(int2str(iYear),'/',month_string,'/',start_day_string);
        end_this_day_search = strcat(int2str(iYear),'/',month_string,'/',end_day_string);

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
                
                disp(strcat('pausing:',datestr(now)));
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

                fileID = fopen(strcat(label,'-',int2str(query_id),'-',int2str(iYear),'-',month_string,'-',int2str(iRetStart),'.xml'),'w','native','utf-8');
                nbytes = fprintf(fileID,'%s',medlineText);
                fclose(fileID);

                % [tree, RootName, DOMnode] = xml_read(strcat(label,'-',int2str(query_id),'-',int2str(iYear),'-',month_string,'-',int2str(iRetStart),'.xml'));
                % save(strcat(label,'-',int2str(query_id),'-',int2str(iYear),'-',month_string,'-',int2str(iRetStart),'.mat'),'tree');

            end
            
        end

    end

end
    
end
   
% hits = regexp(medlineText,'PMID-.*?(?=PMID|</pre>$)','match');
% 
% pmstruct = struct('PubMedID','','PublicationDate','','Title','',...
%                  'Abstract','','Authors','','Citation','');
% 
% % Loop through each article in hits and extract the PubMed ID,
% % publication date, title, abstract, authors, and citation information.
% % Place the information in PMSTRUCT, a MATLAB structure array.
% for n = 1:numel(hits)
%     
%     pmstruct(n).PubMedID          = regexp(hits{n},'(?<=PMID- ).*?(?=\n)','match', 'once');
%     pmstruct(n).PublicationDate   = regexp(hits{n},'(?<=DP  - ).*?(?=\n)','match', 'once');
%     pmstruct(n).Title             = regexp(hits{n},'(?<=TI  - ).*?(?=PG  -|AB  -)','match', 'once');
%     pmstruct(n).Abstract          = regexp(hits{n},'(?<=AB  - ).*?(?=AD  -)','match', 'once');
%     pmstruct(n).Authors           = regexp(hits{n},'(?<=AU  - ).*?(?=\n)','match');
%     pmstruct(n).Citation          = regexp(hits{n},'(?<=SO  - ).*?(?=\n)','match', 'once');
%     
%     pmstruct(n).Language          = regexp(hits{n},'(?<=LA  - ).*?(?=\n)','match', 'once');
%     pmstruct(n).PublicationType   = regexp(hits{n},'(?<=PT  - ).*?(?=\n)','match', 'once');
%     pmstruct(n).DateOfPublication = regexp(hits{n},'(?<=DEP - ).*?(?=\n)','match', 'once');
%     pmstruct(n).TitleAbbreviation = regexp(hits{n},'(?<=TA  - ).*?(?=\n)','match', 'once');
%     pmstruct(n).JournalTitle      = regexp(hits{n},'(?<=JT  - ).*?(?=\n)','match', 'once');
%     pmstruct(n).JournalID         = regexp(hits{n},'(?<=JID - ).*?(?=\n)','match', 'once');
%     pmstruct(n).AID               = regexp(hits{n},'(?<=AID - ).*?(?=\n)','match');
%     
% end


