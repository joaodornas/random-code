
query_label = 'retina';

%%% CONNECT TO DATABASE
conn = database(strcat('references_papers_',query_label),'joaodornas','Senh@1111','Vendor','MySQL','Server','141.44.42.20','PortNumber',7979); 

%%% DEFINE COLUMNS

% getColumnsSelectPubmedArticle = sprintf('SHOW COLUMNS FROM pubmedarticle');

%%% DEFINE QUERIES

% curs = exec(conn,getColumnsSelectPubmedArticle);
% curs = fetch(curs);

%%% PARAMETERS
maxRegisters = 1000;
iStep = 100;

%%% KEYWORDS TABLE
KeywordColumn = 'Keyword';

%%% PUBMEDARTICLE TABLE
TitleColumn = 'ArticleTitle';
% AbstractColumns = {
%     'AbstractEmpty'
%     'AbstractUnassigned'
%     'AbstractUnlabelled'
%     'AbstractComplete'
%     'AbstractIntroduction'
%     'AbstractBackground'
%     'AbstractPurpose'
%     'AbstractObjective'
%     'AbstractAIM'
%     'AbstractMethods'
%     'AbstractSearchMethods'
%     'AbstractSelectionCriteria'
%     'AbstractMaterial'
%     'AbstractMaterialAndMethods'
%     'AbstractDataCollectionAndAnalysis'
%     'AbstractResults'
%     'AbstractMainResults'
%     'AbstractConclusion'
%     'AbstractAuthorsConclusions'
%     'AbstractImportance'
%     'AbstractDesignAndSetting'
%     'AbstractMainOutcomesAndMeasures'
%     'AbstractConclusionsAndRelevance'
%     'AbstractExposure'
%     'AbstractFindings'
%     'AbstractOptions'
%     'AbstractImplications'
%     'AbstractOutcomes'
%     'AbstractInterpretation'
%     'AbstractFunding'
%     'AbstractEvidence'
%     'AbstractContext'
%     'AbstractValues'
%     'AbstractDesign'
%     'AbstractSetting'
%     'AbstractRecommendations'
%     'AbstractMainOutcomeMeasures'
%     'AbstractStatementOfSignificance'
%     'AbstractPatientsAndMethods'
%     'AbstractSubjectsAndMethods'};
AbstractColumn = 'AbstractAllTogether';

%%% GET KEYWORDS FROM KEYWORD TABLE

% disp('...doing Keywords');
% 
% tableName = 'keywords';
% columnName = KeywordColumn;
% 
% getSelectCount = sprintf('SELECT COUNT(*) FROM %s',tableName);
% 
% curs = exec(conn,getSelectCount);
% curs = fetch(curs);
% 
% nTotalRegister = curs.Data{1};
% 
% disp(strcat('nTotalRegister:',int2str(nTotalRegister)));
%             
% if nTotalRegister > iStep
% 
%    nSteps = (nTotalRegister - mod(nTotalRegister,iStep))/iStep;
% 
% else
% 
%    nSteps = 0;
% 
% end
% 
% allSteps = (0:nSteps).*iStep;
% 
% allDataKeyword = [];
% 
% for startStep=allSteps
%     
%     disp(strcat('...',int2str(startStep)));
%     
%     getSelect = sprintf('SELECT %s FROM %s WHERE %s IS NOT NULL LIMIT %s,%s',columnName,tableName,columnName,int2str(startStep),int2str(iStep));
%     
%     curs = exec(conn,getSelect);
%     curs = fetch(curs);
%     
%     allDataKeyword = [allDataKeyword;curs.Data];
%     
% end
% 
% save('Pubmedarticle-keywords.mat','allDataKeyword');

%%% GET KEYWORDS FROM PUBMEDARTICLE TABLE

tableName = 'titleabstract';

getSelectCount = sprintf('SELECT COUNT(*) FROM %s',tableName);

curs = exec(conn,getSelectCount);
curs = fetch(curs);

nTotalRegister = curs.Data{1};
            
if nTotalRegister > iStep

   nSteps = (nTotalRegister - mod(nTotalRegister,iStep))/iStep;

else

   nSteps = 0;

end

allSteps = (0:nSteps).*iStep;

disp('...doing Titles');

allDataTitle = [];

columnName = TitleColumn;

for startStep=allSteps
    
    disp(strcat('...',int2str(startStep)));
    
    getSelect = sprintf('SELECT %s FROM %s WHERE %s IS NOT NULL LIMIT %s,%s',columnName,tableName,columnName,int2str(startStep),int2str(iStep));
    
    curs = exec(conn,getSelect);
    curs = fetch(curs);
    
    allDataTitle = [allDataTitle;curs.Data];
    
end

% [uNouns, uVerbs, uAdverbs, uAdjectives, missing] = getUniqueKeywords(allDataTitle);
% save('Pubmedarticle-title.mat', 'uNouns', 'uVerbs', 'uAdverbs', 'uAdjectives', 'missing');

[uNouns, uVerbs, uAdverbs, uAdjectives, missing] = getUniqueKeywordsNoCleaning(allDataTitle);
save('PubmedarticleNoCleaning-title.mat', 'uNouns', 'uVerbs', 'uAdverbs', 'uAdjectives', 'missing');

disp('...doing Abstracts');


allDataAbstract = [];

columnName = AbstractColumn;

disp(columnName);

for startStep=allSteps

    disp(strcat('...',int2str(startStep)));

    getSelect = sprintf('SELECT %s FROM %s WHERE %s IS NOT NULL LIMIT %s,%s',columnName,tableName,columnName,int2str(startStep),int2str(iStep));

    curs = exec(conn,getSelect);
    curs = fetch(curs);

    allDataAbstract = [allDataAbstract;curs.Data];

end

% [uNouns, uVerbs, uAdverbs, uAdjectives] = getUniqueKeywords(allDataAbstract);
% save('Pubmedarticle-abstract.mat', 'uNouns', 'uVerbs', 'uAdverbs', 'uAdjectives');

[uNouns, uVerbs, uAdverbs, uAdjectives] = getUniqueKeywordsNoCleaning(allDataAbstract);
save('PubmedarticleNoCleaning-abstract.mat', 'uNouns', 'uVerbs', 'uAdverbs', 'uAdjectives');

allData = [];
allData = [allDataTitle, allDataAbstract];

% [uNouns, uVerbs, uAdverbs, uAdjectives] = getUniqueKeywords(allData);
% save('Pubmedarticle-overall.mat', 'uNouns', 'uVerbs', 'uAdverbs', 'uAdjectives');

[uNouns, uVerbs, uAdverbs, uAdjectives] = getUniqueKeywordsNoCleaning(allData);
save('PubmedarticleNoCleaning-overall.mat', 'uNouns', 'uVerbs', 'uAdverbs', 'uAdjectives');

% disp('...doing Abstracts');
% 
% for iAbstract=1:length(AbstractColumn)
%     
%     allDataAbstract(iAbstract).data = [];
% 
%     columnName = AbstractColumn{iAbstract};
%     
%     disp(columnName);
% 
%     for startStep=allSteps
%         
%         disp(strcat('...',int2str(startStep)));
% 
%         getSelect = sprintf('SELECT %s FROM %s WHERE %s IS NOT NULL LIMIT %s,%s',columnName,tableName,columnName,int2str(startStep),int2str(iStep));
% 
%         curs = exec(conn,getSelect);
%         curs = fetch(curs);
% 
%         allDataAbstract(iAbstract).data = [allDataAbstract(iAbstract).data;curs.Data];
% 
%     end
% 
% end
% 
% tmp_data = [];
% for iAbstract=1:length(AbstractColumn)
%     tmp_data = [tmp_data;allDataAbstract(iAbstract).data];
% end
% allDataAbstract = tmp_data;
% 
% [uNouns, uVerbs, uAdverbs, uAdjectives] = getUniqueKeywords(allDataAbstract);
% 
% save('Pubmedarticle-abstract.mat', 'uNouns', 'uVerbs', 'uAdverbs', 'uAdjectives');


