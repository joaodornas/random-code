%%% SEVERAL QUERIES

%%% DATABASE
goal = 'Retina';
version = strcat(int2str(5));
CORE = strcat(goal,'_Healthy_','References');

%%% SERVER
SERVER = '127.0.0.1';
PORT = '8983';
ANALYZER = 'select';
ROWS = strcat('rows=',int2str(0));

%%% DATA (FACET)
%facet.field = 'facet.field=AbstractAllTogether';
facet.field = 'facet.field=PMID';
facet_on = 'facet=on';
indent_on = 'indent=on';
facet.method = 'facet.method=enum';
facet.limit = 'facet.limit=-1';
facet.mincount = 'facet.mincount=1';
facet.sort = 'facet.sort=count';

wt = 'wt=json';
query_ = strcat('http://',SERVER,':',PORT,'/','solr','/',CORE,'/',ANALYZER,'?',ROWS,'&',facet.field,'&',facet_on,'&',indent_on,'&',facet.method,'&',facet.limit,'&',facet.mincount,'&',facet.sort,'&');

options = weboptions('ContentType','json');

%%% OPTION 1: GET RESULTS WITHOUT THE KEYWORDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% q = 'q=';
% 
% d_AbstractAllTogether = '-AbstractAllTogether:';
%     
% load(strcat('__descriptives-',goal,'_v',int2str(version),'.mat'));
% 
% nAnimals = length(animals);
% nPhenomenon = length(phenomenon);
% nNeuron = length(neuron);
% 
% descriptive = [];
% for iAnimals=1:nAnimals
% 
%     descriptive = strcat(descriptive,d_AbstractAllTogether,animals{iAnimals},'%20OR%20');
% 
% end
% for iPhenomenon=1:nPhenomenon
% 
%     descriptive = strcat(descriptive,d_AbstractAllTogether,phenomenon{iPhenomenon},'%20OR%20');
% 
% end
% for iNeuron=1:nNeuron
% 
%     descriptive = strcat(descriptive,d_AbstractAllTogether,neuron{iNeuron},'%20OR%20');
% 
% end
% 
% descriptive = descriptive(1:end-8);
% 
% q = strcat(q,descriptive);
% 
% %%% BUILD QUERY STRING
% 
% query = sprintf('%s%s&%s',query_,q,wt);
% 
% %%% DO QUERY HTTP REQUEST
% 
% whole = urlread(query);
% 
% %%% SAVE RESULT
% 
% fid = fopen('query-lucene-without-keywords.txt','w');
% fprintf(fid,'%s\n',whole);
% fclose(fid);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% OPTION 2: GET RESULTS WITH COMBINATIONS OF KEYWORDS %%%%%%%%%%%%%%%%%%%

q_ = 'q=';

% d_AbstractAllTogether = 'AbstractAllTogether:';
d_AbstractAllTogether = 'ArticleTitle:';
    
load(strcat('__descriptives-',goal,'_v',version,'.mat'));

nAnimals = length(animals);
nPhenomenon = length(phenomenon);
nNeuron = length(neuron);

idx_query = 0;

% only Animals
for iAnimals=1:nAnimals
    
    descriptive = [];

    descriptive = strcat(descriptive,d_AbstractAllTogether,animals{iAnimals});
    
    disp(descriptive);
    
    q = strcat(q_,descriptive);

    %%% BUILD QUERY STRING

    query = sprintf('%s%s&%s',query_,q,wt);

    %%% DO QUERY HTTP REQUEST

    whole = webread(query,options);

    %%% SAVE RESULT
    
    idx_query = idx_query + 1;
    
    all_queries(idx_query).label = strcat('query-lucene-animals-',animals{iAnimals});
    all_queries(idx_query).result = whole;
    
% % %     save(strcat('query-lucene-animals-',animals{iAnimals},'.mat'),'whole');

% % %     fid = fopen(strcat('query-lucene-animals-',animals{iAnimals},'.txt'),'w');
% % %     fprintf(fid,'%s\n',whole);
% % %     fclose(fid);

end

% only Phenomenon
for iPhenomenon=1:nPhenomenon

    descriptive = [];

    descriptive = strcat(descriptive,d_AbstractAllTogether,phenomenon{iPhenomenon});
    
    disp(descriptive);
    
    q = strcat(q_,descriptive);

    %%% BUILD QUERY STRING

    query = sprintf('%s%s&%s',query_,q,wt);

    %%% DO QUERY HTTP REQUEST

    whole = webread(query,options);

    %%% SAVE RESULT
    
    idx_query = idx_query + 1;
    
    all_queries(idx_query).label = strcat('query-lucene-phenomenon-',phenomenon{iPhenomenon});
    all_queries(idx_query).result = whole;
    
% % %     save(strcat('query-lucene-phenomenon-',phenomenon{iPhenomenon},'.mat'),'whole');

% % %     fid = fopen(strcat('query-lucene-phenomenon-',phenomenon{iPhenomenon},'.txt'),'w');
% % %     fprintf(fid,'%s\n',whole);
% % %     fclose(fid);

end

% only Neuron
for iNeuron=1:nNeuron

    descriptive = [];

    descriptive = strcat(descriptive,d_AbstractAllTogether,neuron{iNeuron});
    
    disp(descriptive);
    
    q = strcat(q_,descriptive);

    %%% BUILD QUERY STRING

    query = sprintf('%s%s&%s',query_,q,wt);

    %%% DO QUERY HTTP REQUEST

    whole = webread(query,options);

    %%% SAVE RESULT
    
    idx_query = idx_query + 1;
    
    all_queries(idx_query).label = strcat('query-lucene-neuron-',neuron{iNeuron});
    all_queries(idx_query).result = whole;
    
% % %     save(strcat('query-lucene-neuron-',neuron{iNeuron},'.mat'),'whole');

% % %     fid = fopen(strcat('query-lucene-neuron-',neuron{iNeuron},'.txt'),'w');
% % %     fprintf(fid,'%s\n',whole);
% % %     fclose(fid);

end

% Animal & Neuron
for iAnimals=1:nAnimals
    
    for iNeuron=1:nNeuron
        
        descriptive = [];

        descriptive = strcat(descriptive,d_AbstractAllTogether,animals{iAnimals});  
       
        descriptive = strcat(descriptive,'%20AND%20',d_AbstractAllTogether,neuron{iNeuron});
        
        disp(descriptive);
        
        q = strcat(q_,descriptive);

        %%% BUILD QUERY STRING

        query = sprintf('%s%s&%s',query_,q,wt);

        %%% DO QUERY HTTP REQUEST

        whole = webread(query,options);

        %%% SAVE RESULT
        
        idx_query = idx_query + 1;
    
        all_queries(idx_query).label = strcat('query-lucene-animals-neuron-',animals{iAnimals},'-',neuron{iNeuron});
        all_queries(idx_query).result = whole;  
        
% % % %        save(strcat('query-lucene-animals-neuron-',animals{iAnimals},'-',neuron{iNeuron},'.mat'),'whole');

% % %         fid = fopen(strcat('query-lucene-animals-neuron-',animals{iAnimals},'-',neuron{iNeuron},'.txt'),'w');
% % %         fprintf(fid,'%s\n',whole);
% % %         fclose(fid);
        
    end

end

% Animal & Neuron & Phenomenon
for iAnimals=1:nAnimals
    
    for iNeuron=1:nNeuron
        
        for iPhenomenon=1:nPhenomenon
        
            descriptive = [];

            descriptive = strcat(descriptive,d_AbstractAllTogether,animals{iAnimals});  
       
            descriptive = strcat(descriptive,'%20AND%20',d_AbstractAllTogether,neuron{iNeuron});
            
            descriptive = strcat(descriptive,'%20AND%20',d_AbstractAllTogether,phenomenon{iPhenomenon});
            
            disp(descriptive);
            
            q = strcat(q_,descriptive);

            %%% BUILD QUERY STRING

            query = sprintf('%s%s&%s',query_,q,wt);

            %%% DO QUERY HTTP REQUEST

            whole = webread(query,options);

            %%% SAVE RESULT
            
            idx_query = idx_query + 1;
    
            all_queries(idx_query).label = strcat('query-lucene-animals-neuron-phenomenon-',animals{iAnimals},'-',neuron{iNeuron},'-',phenomenon{iPhenomenon});
            all_queries(idx_query).result = whole;  
            
% % %             save(strcat('query-lucene-animals-neuron-phenomenon-',animals{iAnimals},'-',neuron{iNeuron},'-',phenomenon{iPhenomenon},'.mat'),'whole');

% % %             fid = fopen(strcat('query-lucene-animals-neuron-phenomenon-',animals{iAnimals},'-',neuron{iNeuron},'-',phenomenon{iPhenomenon},'.txt'),'w');
% % %             fprintf(fid,'%s\n',whole);
% % %             fclose(fid);
            
        end
        
    end

end


%%% GET ALL PMIDS from DATABASE

whole = webread(strcat('http://',SERVER,':',PORT,'/','solr','/',CORE,'/',ANALYZER,'?rows=0&indent=on&q=PMID:*&facet=on&facet.field=PMID&facet.method=enum&facet.sort=count&facet.limit=-1&facet.mincount=1&wt=json'),options);

idx_query = idx_query + 1;

all_queries(idx_query).label = strcat('query-lucene-',goal,'-healthy-','PMIDs');
all_queries(idx_query).result = whole;  

% % % save(strcat('query-lucene-',goal,'-healthy-','PMIDs','.mat'),'whole');

% % % fid = fopen(strcat('query-lucene-',goal,'-healthy-','PMIDs','.txt'),'w');
% % % fprintf(fid,'%s\n',whole);
% % % fclose(fid);

save(strcat('query-lucene-',goal,'-v',int2str(version),'.mat'),'all_queries');

