SELECT references_papers_tmp.pubmedarticle.pmid FROM references_papers_tmp.pubmedarticle 
WHERE 
( AbstractEmpty IS NULL 
and  AbstractUnassigned IS NULL
and AbstractUnlabelled IS NULL
and AbstractComplete IS NULL
and AbstractIntroduction IS NULL
and AbstractBackground IS NULL
and AbstractPurpose IS NULL
and AbstractObjective IS NULL
and AbstractAIM IS NULL
and AbstractSearchMethods IS NULL
and AbstractSelectionCriteria IS NULL
and AbstractMaterial IS NULL
and AbstractMaterialAndMethods IS NULL 
and AbstractDataCollectionAndAnalysis IS NULL
and AbstractResults IS NULL
and AbstractMainResults IS NULL
and AbstractConclusion IS NULL
and AbstractAuthorsConclusions IS NULL 
and JournalPubDate > '2016-01-01')