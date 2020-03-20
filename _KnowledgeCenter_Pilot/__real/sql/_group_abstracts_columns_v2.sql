
CREATE TABLE titleabstract_visual (PubmedArticleID int, PMID VARCHAR(45), ArticleTitle TEXT, AbstractAllTogether TEXT, Language VARCHAR(45));

INSERT INTO titleabstract_visual (PubmedArticleID, PMID, ArticleTitle, AbstractAllTogether, Language)
    
SELECT pubmedarticle.idPubmedArticle, pubmedarticle.PMID, pubmedarticle.ArticleTitle,  (SELECT GROUP_CONCAT(DISTINCT AbstractFieldData SEPARATOR ', ') FROM abstractsperpublication WHERE PubmedArticleID = pubmedarticle.idPubmedArticle) , pubmedarticle.Language FROM pubmedarticle

INNER JOIN publicationsperquery ON publicationsperquery.QueryID = 3 AND publicationsperquery.PubmedArticleID = pubmedarticle.idPubmedArticle;
    
# pubmedarticle.Language = 'eng' OR pubmedarticle.Language IS NULL OR pubmedarticle.Language = 'und' AND

