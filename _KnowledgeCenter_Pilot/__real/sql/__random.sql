# SELECT * FROM references_papers_tmp.authors;
# SELECT idAuthors FROM authors WHERE authors.LastName = 'CLAVEL' AND authors.ForeName = '(empty)'
# SELECT idPubmedArticle FROM pubmedarticle WHERE pubmedarticle.PMID = 14084730
# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 20343591
# select * from references_papers_tmp.pubmedarticle order by idPubmedArticle desc limit 10


 ### PUBMED ARTICLES WHERE THE ALGORITHM CRASHED

# SELECT idPubmedArticle FROM pubmedarticle WHERE pubmedarticle.PMID = 3601628
# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 3601628
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 78175 
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 78175

# SELECT idPubmedArticle FROM pubmedarticle WHERE pubmedarticle.PMID = 10576228
# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 10576228
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 78178 
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 78178
 
# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 6243749
 
# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 2752277
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 115984
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 115984
 
# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 10338743

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 16327720
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 115983
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 115983

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 6116519

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 679803
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 125716

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 10528626

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 17306627

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 6322841
# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 7440091
# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 2293070

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 11888395
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 143518
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 143518

 # SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 19829898
 
 ## new
 # SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 6149216
# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 7213232
# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 2251593

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 2487642

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 1496139
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 190675
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 190675

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 15745870
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 190674

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 2926424

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 6086186
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 198950

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 8488080

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 15812663

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 2248284

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 6514501

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 8275902

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 16241019
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 204990
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 204990

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 2174424

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 6499851
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 206103

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 11648255
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 206107

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 15976595
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 206104
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 206104

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 19716260

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 2019739
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 269275
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 269275

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 19011707
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 269277

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 10868188
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 269276
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 269276

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 18766315
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 269414
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 269414

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 9010535
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 269416
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 269416

# SELECT * FROM pubmedarticle WHERE pubmedarticle.PMID = 27738996
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 269747
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 269747

#SHOW variables
#show variables like 'innodb_lock_wait_timeout';

# set net_read_timeout=120;
#SET GLOBAL innodb_lock_wait_timeout = 120;

# SET GLOBAL innodb_buffer_pool_size = 53687091200;
# SET GLOBAL table_open_cache = 80000;
# SET GLOBAL innodb_lock_wait_timeout = 120;

# SET GLOBAL innodb_buffer_pool_dump_now=ON;
 
# SELECT count(*) FROM pubmedarticle WHERE language IS NULL or language = 'eng' or language = 'und'

# SELECT COUNT(*) from titleabstract titabs
#	INNER JOIN pubmedarticle pub
#		ON titabs.PMID = pub.PMID
# WHERE pub.Language IS NULL or pub.Language = 'eng'

# SELECT distinct(language) FROM pubmedarticle
# SELECT COUNT(*) FROM titleabstract
# SELECT *, MATCH (ArticleTitle) AGAINST ('cone rod' IN NATURAL LANGUAGE MODE) AS score FROM titleabstract HAVING score > 0
#Select * from titleabstract limit 0,100
# SELECT * FROM abstractalltogether LIMIT 50000,1000
# SELECT * FROM pubmedarticle WHERE PMID = 12714393
 
# CREATE FULLTEXT INDEX TitleFullText ON titleabstract (ArticleTitle(1024));
# CREATE FULLTEXT INDEX AbstractFullText ON titleabstract (AbstractAllTogether(1024)) ;
# CREATE FULLTEXT INDEX TitleAbstractFullText ON titleabstract (AbstractAllTogether(1024), ArticleTitle(1024));

#UPDATE abstractalltogether aball 
#	INNER JOIN pubmedarticle pub
#		ON aball.PMID = pub.PMID
#SET aball.Title = pub.ArticleTitle
#WHERE pub.PMID IS NOT NULL

#SELECT COUNT(*) FROM wordnet31.senses

# SELECT pubmedarticle.ArticleTitle, MATCH (ArticleTitle) AGAINST ('cone' IN NATURAL LANGUAGE MODE) AS score FROM pubmedarticle HAVING score > 0;
# SELECT *, MATCH (abstractalltogether) AGAINST ('cone' IN NATURAL LANGUAGE MODE) AS score FROM abstractalltogether HAVING score > 0;

#SELECT COUNT(*) FROM titleabstract WHERE ArticleTitle IS NULL

#SELECT COUNT(*) FROM titleabstract WHERE MATCH (ArticleTitle) AGAINST ('retina cone rod photoreceptor bipolar horizontal amacrine ganglion' IN BOOLEAN MODE);
#SELECT COUNT(*) FROM titleabstract WHERE MATCH (AbstractAllTogether) AGAINST ('retina cone rod photoreceptor bipolar horizontal amacrine ganglion' IN BOOLEAN MODE);
#SELECT COUNT(*) FROM titleabstract WHERE MATCH (ArticleTitle,AbstractAllTogether) AGAINST ('retina* cone* rod* photoreceptor* bipolar* horizontal* amacrine* ganglion*' IN BOOLEAN MODE);
  
# SELECT * FROM INFORMATION_SCHEMA.STATISTICS WHERE (table_schema, table_name) = ('references_papers', 'pubmedarticle') ORDER BY seq_in_index;

# SELECT * FROM INFORMATION_SCHEMA.STATISTICS WHERE (table_schema, table_name) = ('references_papers', 'pubmedarticle') AND index_type = 'FULLTEXT' ORDER BY seq_in_index;

# SELECT COUNT(*) FROM pubmedarticle WHERE articleIDPMC IS NOT NULL

# SELECT * FROM pubmedarticle LIMIT 0,1000
# SELECT AbstractComplete, AbstractAllTogether FROM pubmedarticle LIMIT 2000,1000
# SELECT PMID, AbstractAllTogether FROM tmp_abstract WHERE PMID = 8545577

# SELECT * FROM pubmedarticle WHERE PMID = 8545577
# SELECT * FROM pubmedarticle WHERE PMID = 1596219 # 90207
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 90207
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 90207

# SELECT * FROM pubmedarticle WHERE PMID = 15962705
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 134469

# SELECT * FROM pubmedarticle WHERE PMID = 16227880
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 141796

# SELECT * FROM pubmedarticle LIMIT 75040, 2
# SELECT * FROM authorsperpublication WHERE PubmedArticleID = 75041
# SELECT * FROM keywordsperpublication WHERE PubmedArticleID = 75041

 