# SELECT * FROM references_papers_tmp.pubmedarticle WHERE references_papers_tmp.pubmedarticle.JournalPubDate > '2012-01-01';
# SELECT * FROM pubmedarticle WHERE references_papers_tmp.pubmedarticle.JournalVolume >  10

# SELECT  * FROM references_papers_tmp.pubmedarticle WHERE pmid = 2075881

# SELECT  PMID, COUNT(*) FROM pubmedarticle GROUP BY PMID HAVING count(*) = 1;
# SELECT COUNT(*) FROM pubmedarticle WHERE journalpubdate > '2016-01-01';
# SELECT * FROM pubmedarticle WHERE journalpubdate > '2016-01-01';