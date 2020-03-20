
%%% CONNECT TO DATABASE
conn = database('references_papers','joaodornas','Senh@1111','Vendor','MySQL','Server','141.44.42.20','PortNumber',7979); 

AbstractColumns = {
    'AbstractEmpty'
    'AbstractUnassigned'
    'AbstractUnlabelled'
    'AbstractComplete'
    'AbstractIntroduction'
    'AbstractBackground'
    'AbstractPurpose'
    'AbstractObjective'
    'AbstractAIM'
    'AbstractMethods'
    'AbstractSearchMethods'
    'AbstractSelectionCriteria'
    'AbstractMaterial'
    'AbstractMaterialAndMethods'
    'AbstractDataCollectionAndAnalysis'
    'AbstractResults'
    'AbstractMainResults'
    'AbstractConclusion'
    'AbstractAuthorsConclusions'
    'AbstractImportance'
    'AbstractDesignAndSetting'
    'AbstractMainOutcomesAndMeasures'
    'AbstractConclusionsAndRelevance'
    'AbstractExposure'
    'AbstractFindings'
    'AbstractOptions'
    'AbstractImplications'
    'AbstractOutcomes'
    'AbstractInterpretation'
    'AbstractFunding'
    'AbstractEvidence'
    'AbstractContext'
    'AbstractValues'
    'AbstractDesign'
    'AbstractSetting'
    'AbstractRecommendations'
    'AbstractMainOutcomeMeasures'
    'AbstractStatementOfSignificance'
    'AbstractPatientsAndMethods'
    'AbstractSubjectsAndMethods'};

nAbstractColumns = length(AbstractColumns);

select = 'SELECT PMID, ';
select_abstract = select;

for iAbstract=1:nAbstractColumns
    
    select_abstract = sprintf('%s%s, ',select_abstract,AbstractColumns{iAbstract});
    
end

idx_comma = strfind(select_abstract,',');
select_abstract(idx_comma(end)) = [];

select_count = 'SELECT COUNT(*) FROM pubmedarticle';
curs = exec(conn,select_count);
curs = fetch(curs);

nTotalPapers = curs.Data{1};
nSteps = (nTotalPapers - mod(nTotalPapers,1000))/1000;
Steps = (0:nSteps+1).*1000;

for iStep=1:length(Steps)
    
    disp(int2str(iStep));

    select = sprintf('%s FROM pubmedarticle LIMIT %s,%s',select_abstract,int2str(Steps(iStep)),int2str(1000));

    setdbprefs('NullStringWrite','<null>');
    setdbprefs('NullStringRead','<null>');

    curs = exec(conn,select);
    curs = fetch(curs);
    
    allData = curs.Data;
    
    nPapers = size(allData,1);
    
    for iPaper=1:nPapers
        
        getData = allData(iPaper,2:end);
        PMID = allData(iPaper,1);
        
        getData = strjoin(getData);
        getDatacell = {getData};
        
        WHERECLAUSE = {sprintf('where PMID = %s',int2str(PMID{1}))};
        COLUMNS = {'AbstractAllTogether'};
        
        update(conn,'pubmedarticle',COLUMNS,getDatacell,WHERECLAUSE);
        
        clear getData PMID
        
    end
       
    close(curs);
    
    close(conn);
    
    conn = database('references_papers','joaodornas','Senh@1111','Vendor','MySQL','Server','141.44.42.20','PortNumber',7979);
    
    clear allData 
    
end

close(conn);
    
    