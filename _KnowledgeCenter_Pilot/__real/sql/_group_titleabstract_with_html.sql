
# CREATE TABLE titleabstract (PMID VARCHAR(45), ArticleTitle TEXT, AbstractAllTogether TEXT);

INSERT INTO titleabstract(PMID, ArticleTitle, AbstractAllTogether)

SELECT PMID, ArticleTitle, AbstractAllTogether from references_papers_retina.titleabstract 
    
#SET SQL_SAFE_UPDATES = 0;
#set innodb_lock_wait_timeout=10000;
#set net_read_timeout=120;

#UPDATE pubmedarticle pub 
#	INNER JOIN tmp_abstract tmp
#		ON pub.PMID = tmp.PMID
#SET pub.AbstractAllTogether = tmp.AbstractAllTogether
#WHERE pub.idPubmedArticle IS NOT NULL
