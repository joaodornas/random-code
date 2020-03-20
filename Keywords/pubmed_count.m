% pubmed_trend.m
term=input('Term:','s');
for j=1990:2011
query=['http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&amp;term=' urlencode(term) '[tiab]%20AND%20' num2str(j) '[dp]&amp;rettype=count'];
 
docNode = xmlread(query);
 
count(j)=str2double(docNode.getElementsByTagName('Count').item(0).getFirstChild.getNodeValue);
 
end;
figure('name',['occurences of the term ' term ' in Pubmed (in title or abstract'])
plot(count(1990:2011));
set(gca,'xtick', 1:5:45,'xtickl',1990:5:2010);
title(term);
ylabel 'number of articles'
xlabel 'year'