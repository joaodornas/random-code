function textData = getUniqueKeywordsFromLucene(Field,CORE)

%%% SERVER
SERVER = '35.190.184.245';
PORT = '8983';
ANALYZER = 'select';
ROWS = strcat('rows=',int2str(0));

%%% DATA (FACET)
facet.field = strcat('facet.field=',Field);
facet_on = 'facet=on';
indent_on = 'indent=on';
facet.method = 'facet.method=enum';
facet.limit = 'facet.limit=-1';
facet.mincount = 'facet.mincount=1';
facet.sort = 'facet.sort=count';

query_ = strcat('http://',SERVER,':',PORT,'/','solr','/',CORE,'/',ANALYZER,'?',ROWS,'&',facet.field,'&',facet_on,'&',indent_on,'&',facet.method,'&',facet.limit,'&',facet.mincount,'&',facet.sort,'&');
q_ = strcat('q=',Field,':*');
wt = 'wt=json';

query = sprintf('%s%s&%s',query_,q_,wt);

options = weboptions('ContentType','json');
textData = webread(query,options);


end