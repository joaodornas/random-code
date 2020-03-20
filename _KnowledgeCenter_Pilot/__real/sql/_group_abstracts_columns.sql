
CREATE TABLE titleabstract (PubmedArticleID int, PMID VARCHAR(45), ArticleTitle TEXT, AbstractAllTogether TEXT);

INSERT INTO titleabstract (PubmedArticleID, PMID, ArticleTitle, AbstractAllTogether)

SELECT idPubmedArticle, PMID, ArticleTitle,
	CONCAT_WS(', ', AbstractEmpty,
   AbstractUnassigned,
   AbstractUnlabelled,
   AbstractComplete,
   AbstractIntroduction,
    AbstractBackground,
    AbstractPurpose,
    AbstractObjective,
    AbstractAIM,
    AbstractMethods,
    AbstractSearchMethods,
    AbstractSelectionCriteria,
    AbstractMaterial,
    AbstractMaterialAndMethods,
    AbstractDataCollectionAndAnalysis,
    AbstractResults,
    AbstractMainResults,
    AbstractConclusion,
    AbstractAuthorsConclusions,
    AbstractImportance,
    AbstractDesignAndSetting,
    AbstractMainOutcomesAndMeasures,
    AbstractConclusionsAndRelevance,
    AbstractExposure,
    AbstractFindings,
    AbstractOptions,
    AbstractImplications,
    AbstractOutcomes,
    AbstractInterpretation,
    AbstractFunding,
    AbstractEvidence,
    AbstractContext,
    AbstractValues,
    AbstractDesign,
    AbstractSetting,
    AbstractRecommendations,
    AbstractMainOutcomeMeasures,
    AbstractStatementOfSignificance,
    AbstractPatientsAndMethods,
    AbstractSubjectsAndMethods,
    AbstractBackgroundAndMethods,
    AbstractMethodsAndResults,
    AbstractBenefitsCostsAndHarms,
    AbstractValidation,
    AbstractSponsor) from pubmedarticle WHERE pubmedarticle.Language = 'eng' OR pubmedarticle.Language IS NULL OR pubmedarticle.Language = 'und'
    
#SET SQL_SAFE_UPDATES = 0;
#set innodb_lock_wait_timeout=10000;
#set net_read_timeout=120;

#UPDATE pubmedarticle pub 
#	INNER JOIN tmp_abstract tmp
#		ON pub.PMID = tmp.PMID
#SET pub.AbstractAllTogether = tmp.AbstractAllTogether
#WHERE pub.idPubmedArticle IS NOT NULL
