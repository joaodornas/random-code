
function my_pubmed2mysql(which_query,start_year,end_year)

switch which_query
    
    case 1
        
        query_label = 'retina';    

        query_id = 2;

        I_Am_Testing = 0;
        
        OS = 'MAC';

        pubmed2mysql(query_label,query_id,start_year,end_year,I_Am_Testing,OS);
        
end

function pubmed2mysql(query_label,query_id,start_year,end_year,I_Am_Testing,OS)
    

switch OS
    
     
    case 'MAC'

        full_path = '/Users/dornas/Download/RETINA';
       
        log_folder = '/Users/dornas/Download/_LOGS/';
        
        javaaddpath('/Users/joaodornas/Dropbox (joaodornas)/MindScience/_CODE/reference-database/Reference-Database/__tools/DB_drivers/mysql-connector-java-5.1.40-bin.jar');
        
    case 'GCLOUD'
        
        full_path = '/home/joaodornas/_PUBMED-references';
        
        log_folder = '/home/joaodornas/_LOGS/';
        
        javaaddpath('/home/joaodornas/DB_drivers/mysql-connector-java-5.1.40-bin.jar');
        
end

full_path_year = strcat(full_path,'/',query_label,'-',int2str(query_id),'/');

list = dir(strcat(full_path_year,'*.xml'));

% for iYear=start_year:end_year
      
    if ~isempty(list)
        
       nFiles = length(list);
       
       for iFile=1:nFiles
           
           name = list(iFile).name;
           idx = strfind(name,'-');
          
           iYear = str2num(list(iFile).name(idx(2)+1:idx(2)+4));
       
       fileID = fopen(strcat(log_folder,'log','-',query_label,'-',int2str(query_id),'-',int2str(start_year),'-',int2str(end_year),'.txt'),'a+');
       
       if iYear >= start_year 
           
       try
           
          [xml_tree, RootName, DOMnode] = xml_read(strcat(full_path_year,list(iFile).name));
          
          fprintf(fileID,'%s\n',strcat('FileLoaded:',int2str(iYear),'-',int2str(iFile)));
          
       catch ME
           
           msg = getReport(ME);
           
           disp(strcat('ErrorLoadingFile:',int2str(iYear),'-',int2str(iFile)));
           
           disp(msg);
           
           fprintf(fileID,'%s\n',strcat('ErrorLoadingFile:',int2str(iYear),'-',int2str(iFile)));
           
           fprintf(fileID,'%s\n',msg);
           
       end
          
           if ~isempty(xml_tree) && isstruct(xml_tree)
              
              apriori_fields = fieldnames(xml_tree);
              
              if ~strcmp('ERROR',apriori_fields)
          
                  nPublications = length(xml_tree.PubmedArticle);

                  if nPublications > 0
                      
                      if I_Am_Testing
                        conn = database('references_papers_tmp','joaodornas','senha1111','Vendor','MySQL','Server','127.0.0.1'); %%% TEST DATABASE
                      else
                        conn = database(strcat('references_papers_',query_label),'root','senha1111','Vendor','MySQL','Server','192.168.0.22','PortNumber',3306);
                      end

                      set(conn,'AutoCommit','off');

                      for iPub=1:nPublications
                          
                          disp(strcat(int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));
                          
                          clearvars -except iPub iYear nPublications xml_tree apriori_fields iFile nFiles list full_path_year start_year end_year full_path conn query_id query_label iPub_data I_Am_Testing fileID OS log_folder
                          
                          %%%%% LOAD DATA ABOUT PUBLICATION FROM STRUCTURE
                          
                          %%% MedlineCitation
                          
                          fields_MedlineCitation = fieldnames(xml_tree.PubmedArticle(iPub).MedlineCitation);
                          
                          if sum(strcmp('PMID',fields_MedlineCitation)) ~= 0
                          
                                PMID = xml_tree.PubmedArticle(iPub).MedlineCitation.PMID.CONTENT;
                               
                          else
                              
                                PMID = cell.empty;
                              
                          end
                          
                          if sum(strcmp('DateCreated',fields_MedlineCitation)) ~= 0

                                Year = xml_tree.PubmedArticle(iPub).MedlineCitation.DateCreated.Year;
                                Month = xml_tree.PubmedArticle(iPub).MedlineCitation.DateCreated.Month;
                                Day = xml_tree.PubmedArticle(iPub).MedlineCitation.DateCreated.Day;
                                
                                DateCreated = {datestr(datenum(strcat(Day,'.',Month,'.',Year),'dd.mm.yyyy'),26)};
                                
                          else
                              
                                DateCreated = cell.empty;
                              
%                               DateCreated.Year = cell.empty;
%                               DateCreated.Month = cell.empty;
%                               DateCreated.Day = cell.empty;
                              
                          end
                          
                          if sum(strcmp('DateRevised',fields_MedlineCitation)) ~= 0

                                Year = xml_tree.PubmedArticle(iPub).MedlineCitation.DateRevised.Year;
                                Month = xml_tree.PubmedArticle(iPub).MedlineCitation.DateRevised.Month;
                                Day = xml_tree.PubmedArticle(iPub).MedlineCitation.DateRevised.Day;
                                
                                DateRevised = {datestr(datenum(strcat(Day,'.',Month,'.',Year),'dd.mm.yyyy'),26)};
                                
                          else
                              
                                DateRevised = cell.empty;
                              
%                               DateRevised.Year = cell.empty;
%                               DateRevised.Month = cell.empty;
%                               DateRevised.Day = cell.empty;
                              
                          end
                          
                          if sum(strcmp('DateCompleted',fields_MedlineCitation)) ~= 0

                                Year = xml_tree.PubmedArticle(iPub).MedlineCitation.DateCompleted.Year;
                                Month = xml_tree.PubmedArticle(iPub).MedlineCitation.DateCompleted.Month;
                                Day = xml_tree.PubmedArticle(iPub).MedlineCitation.DateCompleted.Day;
                                
                                DateCompleted = {datestr(datenum(strcat(Day,'.',Month,'.',Year),'dd.mm.yyyy'),26)};
                                
                          else
                              
                              DateCompleted = cell.empty;
                              
%                               DateCompleted.Year = cell.empty;
%                               DateCompleted.Month = cell.empty;
%                               DateCompleted.Day = cell.empty;
                              
                          end  
                        
                          %%% MedlineCitation - Article
                          
                          fields_Article = fieldnames(xml_tree.PubmedArticle(iPub).MedlineCitation.Article);
                          
                          if sum(strcmp('PublicationTypeList',fields_Article)) ~= 0
                          
                                PublicationType = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.PublicationTypeList.PublicationType.CONTENT;
                          
                          else
                              
                              PublicationType = cell.empty;
                              
                          end
                          
                          if sum(strcmp('ArticleTitle',fields_Article)) ~= 0
                              
                                ArticleTitle = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.ArticleTitle;
                                
                          else
                              
                              ArticleTitle = cell.empty;
                              
                          end
                          
                          if sum(strcmp('ArticleDate',fields_Article)) ~= 0

                                Year = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.ArticleDate.Year;
                                Month = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.ArticleDate.Month;
                                Day = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.ArticleDate.Day;
                                
                                ArticleDate = {datestr(datenum(strcat(Day,'.',Month,'.',Year),'dd.mm.yyyy'),26)};
                                
                          else
                              
                              ArticleDate = cell.empty;
                              
%                               ArticleDate.Year = cell.empty;
%                               ArticleDate.Month = cell.empty;
%                               ArticleDate.Day = cell.empty;
                              
                          end
                          
                          if sum(strcmp('Pagination',fields_Article)) ~= 0
                              
                                Pagination = num2str(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Pagination.MedlinePgn);
                                
                          else
                              
                              Pagination = cell.empty;
                              
                          end
                          
                          if sum(strcmp('Language',fields_Article)) ~= 0
                          
                                Language = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Language;
                                
                          else
                              
                              Language = cell.empty;
                              
                          end
                          
                          %%% PubmedData - ArticleIdList
                          
                          fields_ArticleIdList = fieldnames(xml_tree.PubmedArticle(iPub).PubmedData.ArticleIdList);
                          
                          articleID_pubmed = cell.empty;
                          articleID_doi = cell.empty;
                          articleID_pii = cell.empty;
                          articleID_pmc = cell.empty;
                          articleID_mid = cell.empty;
                          articleID_pmcid = cell.empty;
                          
                          if sum(strcmp('ArticleId',fields_ArticleIdList)) ~= 0
                              
                              for iArt=1:length(xml_tree.PubmedArticle(iPub).PubmedData.ArticleIdList.ArticleId)
                                  
                                IdType = xml_tree.PubmedArticle(iPub).PubmedData.ArticleIdList.ArticleId(iArt).ATTRIBUTE.IdType;
                                
                                if iscell(IdType)

                                    eval(strcat('articleID_',IdType,'{1} = xml_tree.PubmedArticle(iPub).PubmedData.ArticleIdList.ArticleId(',int2str(iArt),').CONTENT;'));
                                
                                    eval(strcat('articleID_',IdType,'{1} = num2str(articleID_',IdType,'{1});'));
                                    
                                else
                                    
                                    eval(strcat('articleID_',IdType,' = xml_tree.PubmedArticle(iPub).PubmedData.ArticleIdList.ArticleId(',int2str(iArt),').CONTENT;'));
                                
                                    eval(strcat('articleID_',IdType,' = num2str(articleID_',IdType,');'));
                                    
                                end
                                
                              end
                 
                          end
                          
                          %%% MedlineCitation - Article - Journal
                          
                          fields_Journal = fieldnames(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Journal);
                          
                          if sum(strcmp('Title',fields_Journal)) ~= 0

                                JournalTitle = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Journal.Title;
                                
                                JournalTitle = strrep(JournalTitle,'"','');
                                JournalTitle = strrep(JournalTitle,'''','');
                             
                          else
                              
                              JournalTitle = cell.empty;
                              
                          end
                          
                          if sum(strcmp('ISOAbbreviation',fields_Journal)) ~= 0
                          
                                JournalISOAbbreviation = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Journal.ISOAbbreviation;
                                
                          else
                              
                              JournalISOAbbreviation = cell.empty;
                              
                          end
                          
                          if sum(strcmp('JournalIssue',fields_Journal)) ~= 0
                              
                              fields_JournalIssue = fields(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Journal.JournalIssue);
                              
                              if sum(strcmp('Volume',fields_JournalIssue)) ~= 0

                                    JournalVolume = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Journal.JournalIssue.Volume;
                                
                              else
                                  
                                    JournalVolume = cell.empty;
                                  
                              end
                              
                              if sum(strcmp('Issue',fields_JournalIssue)) ~= 0

                                    JournalIssue = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Journal.JournalIssue.Issue;
                                
                              else
                                  
                                    JournalIssue = cell.empty;
                                  
                              end
                              
                              if sum(strcmp('PubDate',fields_JournalIssue)) ~= 0
                                
                                    fields_PubDate = fieldnames(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Journal.JournalIssue.PubDate);
                                
                                    if sum(strcmp('Year',fields_PubDate)) ~= 0; Year = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Journal.JournalIssue.PubDate.Year; else Year = '1'; end
                                
                                    Day = '1';
                                    Month = '1';
                                    
                                    JournalPubDate = {datestr(datenum(strcat(Day,'.',Month,'.',Year),'dd.mm.yyyy'),26)};
                                    
                              else
                                  
                                  JournalPubDate = cell.empty;
                                  
                              end
                                
                          else
                              
                              JournalVolume = cell.empty;
                              JournalPubDate = cell.empty;
                              JournalIssue = cell.empty;
                              
                          end
                          
                          %%% MedlineCitation - MedlineJournalInfo
                          
                          fields_MedlineJournalInfo = fieldnames(xml_tree.PubmedArticle(iPub).MedlineCitation.MedlineJournalInfo);

                          if sum(strcmp('Country',fields_MedlineJournalInfo)) ~= 0
                              
                                JournalCountry = xml_tree.PubmedArticle(iPub).MedlineCitation.MedlineJournalInfo.Country;
                          
                          else
                              
                              JournalCountry = cell.empty;
                              
                          end
                          
                          if sum(strcmp('MedlineTA',fields_MedlineJournalInfo)) ~= 0
                              
                                JournalMedlineTA = xml_tree.PubmedArticle(iPub).MedlineCitation.MedlineJournalInfo.MedlineTA;
                                
                          else
                              
                              JournalMedlineTA = cell.empty;
                              
                          end
                          
                          if sum(strcmp('NlmUniqueID',fields_MedlineJournalInfo)) ~= 0
                              
                                JournalNlmUniqueID = xml_tree.PubmedArticle(iPub).MedlineCitation.MedlineJournalInfo.NlmUniqueID;
                                
                                if isnumeric(JournalNlmUniqueID); JournalNlmUniqueID = num2str(JournalNlmUniqueID); else JournalNlmUniqueID = char(JournalNlmUniqueID); end
                                
                          else
                              
                              JournalNlmUniqueID = cell.empty;
                              
                          end
                          
                          if sum(strcmp('ISSNLinking',fields_MedlineJournalInfo)) ~= 0
                              
                                JournalISSNLinking = xml_tree.PubmedArticle(iPub).MedlineCitation.MedlineJournalInfo.ISSNLinking;
                                
                                if isnumeric(JournalISSNLinking); JournalISSNLinking = num2str(JournalISSNLinking); else JournalISSNLinking = char(JournalISSNLinking); end
                                
                          else
                              
                              JournalISSNLinking = cell.empty;
                              
                          end

                          % idJournal =

                          %%% MedlineCitation - Article - Abstract
                          
                          Abstract_empty = cell.empty;
                          Abstract_unassigned = cell.empty;
                          Abstract_unlabelled = cell.empty;
                          Abstract_complete = cell.empty;
                          Abstract_introduction = cell.empty;
                          Abstract_background = cell.empty;
                          Abstract_purpose = cell.empty;
                          Abstract_objective = cell.empty;
                          Abstract_aim = cell.empty;
                          Abstract_methods = cell.empty;
                          Abstract_search_methods = cell.empty;
                          Abstract_selection_criteria = cell.empty;
                          Abstract_material = cell.empty;
                          Abstract_material_and_methods = cell.empty;
                          Abstract_data_collection_and_analysis = cell.empty;
                          Abstract_results = cell.empty;
                          Abstract_main_results = cell.empty;
                          Abstract_conclusion = cell.empty;
                          Abstract_authors_conclusions = cell.empty;
                          Abstract_importance = cell.empty;
                          Abstract_design_and_setting = cell.empty;
                          Abstract_exposure = cell.empty;
                          Abstract_main_outcomes_and_measures = cell.empty;
                          Abstract_conclusions_and_relevance = cell.empty;
                          Abstract_findings = cell.empty;
                          Abstract_options = cell.empty;
                          Abstract_implications = cell.empty;
                          Abstract_outcomes = cell.empty;
                          Abstract_interpretation = cell.empty;
                          Abstract_funding = cell.empty;
                          Abstract_evidence = cell.empty;
                          Abstract_context = cell.empty;
                          Abstract_values = cell.empty;
                          Abstract_design = cell.empty;
                          Abstract_setting = cell.empty;
                          Abstract_recommendations = cell.empty;
                          Abstract_main_outcome_measures = cell.empty;
                          Abstract_statement_of_significance = cell.empty;
                          Abstract_patients_and_methods = cell.empty;
                          Abstract_subjects_and_methods = cell.empty;
                          Abstract_background_and_methods = cell.empty;
                          Abstract_methods_and_results = cell.empty;
                          Abstract_benefits_costs_and_harms = cell.empty;
                          Abstract_validation = cell.empty;
                          Abstract_sponsor = cell.empty;
                          
                          if (sum(strcmp('Abstract',fields_Article)) ~= 0)
                          
                              fields_Abstract = fieldnames(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract); 

                              if (sum(strcmp('AbstractText',fields_Abstract)) ~= 0) && (isstruct(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText)) && (length(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText) > 1)

                                    Abstract_complete = cell.empty;

                                    for iAb=1:length(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText)
                                        
                                        fields_ATTRIBUTE = fieldnames(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText(iAb).ATTRIBUTE);
                                        
                                        if (sum(strcmp('NlmCategory',fields_ATTRIBUTE)) ~= 0)

                                            category = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText(iAb).ATTRIBUTE.NlmCategory;
                                            
                                        else
                                            
                                            if (sum(strcmp('Label',fields_ATTRIBUTE)) ~= 0)
                                            
                                                category = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText(iAb).ATTRIBUTE.Label;
                                            
                                            end
                                            
                                        end
                                        
                                        if ~isempty(category)
                                        
                                            category = strrep(category,' ','_');
                                            category = strrep(category,'''','');
    
                                            if strcmp('conclusions',lower(category)); category = 'conclusion'; end
                                            if strcmp('objectives',lower(category)); category = 'objective'; end
                                            if strcmp('method',lower(category)); category = 'methods'; end  
                                          
                                            eval(strcat('Abstract_',lower(category),' = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText(',int2str(iAb),').CONTENT;'));

                                        else
                                        
                                            %%% 'Empty Caterory'.
                                            
                                            category = 'empty';
                                            
                                            eval(strcat('Abstract_',lower(category),' = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText(',int2str(iAb),').CONTENT;'));
                                            
                                        end
                                        
                                    end

                              elseif length(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText) == 1

                                  if isstruct(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText)
                                      
                                        Abstract_complete = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText.CONTENT;

                                  else
                                  
                                      Abstract_complete = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText;
                                      
                                  end
                                  
                              elseif ischar(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText)
                                  
                                  Abstract_complete = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract.AbstractText;
                                  
                              end

                          end
                          
                          WHOLE_Abstract = strcat( ' <empty> ', Abstract_empty, ...
                           ' <unassigned> ', Abstract_unassigned, ...
                           ' <unlabelled> ', Abstract_unlabelled, ...
                           ' <complete> ', Abstract_complete, ...
                           ' <introduction> ', Abstract_introduction, ...
                           ' <background> ', Abstract_background, ...
                           ' <purpose> ', Abstract_purpose, ...
                           ' <objective> ', Abstract_objective, ...
                           ' <aim> ', Abstract_aim, ...
                           ' <methods> ', Abstract_methods, ...
                           ' <search_methods> ', Abstract_search_methods, ...
                           ' <selection_criteria> ', Abstract_selection_criteria, ...
                           ' <material> ', Abstract_material, ...
                           ' <material_and_methods> ', Abstract_material_and_methods, ...
                           ' <data_collection_and_analysis> ', Abstract_data_collection_and_analysis, ...
                           ' <results> ', Abstract_results, ...
                           ' <main_results> ', Abstract_main_results, ...
                           ' <conclusion> ', Abstract_conclusion, ...
                           ' <authors_conclusions> ', Abstract_authors_conclusions, ...
                           ' <importance> ', Abstract_importance, ...
                           ' <design_and_setting> ', Abstract_design_and_setting, ...
                           ' <exposure> ', Abstract_exposure, ...
                           ' <main_outcomes_and_measures> ', Abstract_main_outcomes_and_measures, ...
                           ' <conclusions_and_relevance> ', Abstract_conclusions_and_relevance, ...
                           ' <findings> ', Abstract_findings, ...
                           ' <options> ', Abstract_options, ...
                           ' <implications> ', Abstract_implications, ...
                           ' <outcomes> ', Abstract_outcomes, ...
                           ' <interpretation> ', Abstract_interpretation, ...
                           ' <funding> ', Abstract_funding, ...
                           ' <evidence> ', Abstract_evidence, ...
                           ' <context> ', Abstract_context, ...
                           ' <values> ', Abstract_values, ...
                           ' <design> ', Abstract_design, ...
                           ' <setting> ', Abstract_setting, ...
                           ' <recommendations> ', Abstract_recommendations, ...
                           ' <main_outcome_measures> ', Abstract_main_outcome_measures, ...
                           ' <statement_of_significance> ', Abstract_statement_of_significance, ...
                           ' <patients_and_methods> ', Abstract_patients_and_methods, ...
                           ' <subjects_and_methods> ', Abstract_subjects_and_methods, ...
                           ' <background_and_methods> ', Abstract_background_and_methods, ...
                           ' <methods_and_results> ', Abstract_methods_and_results, ...
                           ' <benefits_costs_and_harms> ', Abstract_benefits_costs_and_harms, ...
                           ' <validation> ', Abstract_validation, ...
                           ' <sponsor> ', Abstract_sponsor ...
                          );
                          
                          %%% MedlineCitation - Article - AuthorList
                          
                          if (sum(strcmp('AuthorList',fields_Article)) ~= 0)
                          
                                fields_AuthorList = fieldnames(xml_tree.PubmedArticle(iPub).MedlineCitation.Article.AuthorList); 
                          
                                if sum(strcmp('Author',fields_AuthorList)) ~= 0

                                    All_Authors = xml_tree.PubmedArticle(iPub).MedlineCitation.Article.AuthorList.Author;
                                    
                                    nAuthors = length(All_Authors);
                                    
                                    for iA=1:nAuthors
                                        
                                        fields_Authors = fieldnames(All_Authors(iA));
                                        
                                        if (sum(strcmp('LastName',fields_Authors)) == 0); All_Authors(iA).LastName = 'none'; end
                                        if (sum(strcmp('ForeName',fields_Authors)) == 0); All_Authors(iA).ForeName = 'none'; end
                                        if (sum(strcmp('Initials',fields_Authors)) == 0); All_Authors(iA).Initials = 'none'; end
                                        
                                    end
                                
                                else
                              
                                    All_Authors = cell.empty;
                              
                                end
                                
                          else
                              
                              All_Authors = cell.empty;
                              
                          end
                          
                          %%% MedlineCitation - KeywordList
                          
                          if sum(strcmp('KeywordList',fields_MedlineCitation)) ~= 0
                          
                                fields_KeywordList = fieldnames(xml_tree.PubmedArticle(iPub).MedlineCitation.KeywordList);
                          
                                if sum(strcmp('Keyword',fields_KeywordList)) ~= 0

                                    All_Keywords = xml_tree.PubmedArticle(iPub).MedlineCitation.KeywordList.Keyword;
                                
                                else
                              
                                All_Keywords = cell.empty;
                              
                                end
                                
                          else
                              
                              All_Keywords = cell.empty;
                              
                          end

                          %%%%% LOAD DATA ABOUT PUBLICATION TO DATABASE
                          
                          if ~isempty(PMID)
                             
                                %% conn = checkConn(conn);
                                
                                select_check_PMID = strcat('SELECT COUNT(*) FROM pubmedarticle WHERE pubmedarticle.PMID = ',PMID);
                                
% %                                 curs = exec(conn,select_check_PMID);
% %                                 curs = fetch(curs);

                                curs = my_exec(conn,select_check_PMID);
        
                                if curs.Data{1} == 0; PMID_exist = 0; else PMID_exist = 1; end

                                close(curs);
                                
                                if ~PMID_exist
                                    
                                    %%% CHECK IF JOURNAL RECORD EXIST
                                    
                                    if ~isempty(JournalTitle)
                                       
                                        %% conn = checkConn(conn);
                                        
                                        select_check_JournalTitle = strcat('SELECT COUNT(*) FROM journal WHERE journal.Title = "',JournalTitle,'"');
                                        
%                                         curs = exec(conn,select_check_JournalTitle);
%                                         curs = fetch(curs);

                                        curs = my_exec(conn,select_check_JournalTitle);

                                        if curs.Data{1} == 0; JournalTitle_exist = 0; else JournalTitle_exist = 1; end
                                
                                        close(curs);
                                        
                                        if JournalTitle_exist
                                            
                                            %% conn = checkConn(conn);
                                            
                                            select_check_JournalTitle = strcat('SELECT idJournal FROM journal WHERE journal.Title = "',JournalTitle,'"');
                                        
%                                             curs = exec(conn,select_check_JournalTitle);
%                                             curs = fetch(curs);

                                            curs = my_exec(conn,select_check_JournalTitle);

                                            idJournal = curs.Data{1};
                                        
                                            close(curs);
                                            
                                        else
                                           
                                            table_name = 'journal';
                                            
                                            data{1} = JournalTitle; colnames{1} = 'Title';
                                            if ~isempty(JournalISOAbbreviation); data{end+1} = JournalISOAbbreviation; colnames{end+1} = 'ISOAbbreviation'; end
                                            if ~isempty(JournalMedlineTA); data{end+1} = JournalMedlineTA; colnames{end+1} = 'MedlineTA'; end
                                            if ~isempty(JournalNlmUniqueID); data{end+1} = JournalNlmUniqueID; colnames{end+1} = 'NlmUniqueID'; end
                                            if ~isempty(JournalISSNLinking); data{end+1} = JournalISSNLinking; colnames{end+1} = 'ISSNLinking'; end
                                            if ~isempty(JournalCountry); data{end+1} = JournalCountry; colnames{end+1} = 'Country'; end
                                            
                                            try
                                                
                                                %% conn = checkConn(conn);
                                                
                                                fastinsert(conn,table_name,colnames,data);
                                                clear colnames data
                                                
                                                my_commit(conn);
                            
                                                disp(strcat('Commit:',int2str(iYear),'-',int2str(iFile)));
                            
                                                fprintf(fileID,'%s\n',strcat('Commit:',table_name,':',int2str(iYear),'-',int2str(iFile)));
               
                                                
                                            catch ME
                                                
                                                rollback(conn);
                                                
                                                msg = getReport(ME);
                                                
                                                disp(msg);
                                                
                                                disp(strcat('ErrorOnFastInsert:',table_name,':',int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));
                                                
                                                disp(msg);
                                                
                                                fprintf(fileID,'%s\n',strcat('ErrorOnFastInsert:',table_name,':',int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));
                                               
                                                fprintf(fileID,'%s\n',msg);
                                                
                                                disp(strcat('Rollback:',int2str(iYear),'-',int2str(iFile)));
                            
                                                fprintf(fileID,'%s\n',strcat('Rollback:',int2str(iYear),'-',int2str(iFile)));
                                                
                                                break
                                                
                                            end
                                            
                                            %% conn = checkConn(conn);
                                            
                                            select_check_JournalTitle = strcat('SELECT idJournal FROM journal WHERE journal.Title = "',JournalTitle,'"');
                                         
%                                             curs = exec(conn,select_check_JournalTitle);
%                                             curs = fetch(curs);

                                            curs = my_exec(conn,select_check_JournalTitle);

                                            idJournal = curs.Data{1};
                                            
                                            close(curs);
                                            
                                        end
                                        
                                    else
                                        
                                        idJournal = cell.empty;
                                        
                                    end
                                    
                                    %%% INSERT ARTICLE DATA
                                    table_name = 'pubmedarticle';
                                    
                                    data{1} = str2num(PMID);
                                    colnames{1} = 'PMID';
                                    
                                    if ~isempty(DateCreated); data{end+1} = DateCreated; colnames{end+1} = 'DateCreated'; end
                                    
                                    if ~isempty(DateRevised); data{end+1} = DateRevised; colnames{end+1} = 'DateRevised'; end
                                    
                                    if ~isempty(DateCompleted); data{end+1} = DateCompleted; colnames{end+1} = 'DateCompleted'; end
                                    
                                    data{end+1} = idJournal;
                                    colnames{end+1} = 'idJournal';
                                    
                                    if ~isempty(JournalVolume); data{end+1} = JournalVolume; colnames{end+1} = 'JournalVolume'; end
                                    
                                    if ~isempty(JournalIssue); data{end+1} = JournalIssue; colnames{end+1} = 'JournalIssue'; end
                                    
                                    if ~isempty(JournalPubDate); data{end+1} = JournalPubDate; colnames{end+1} = 'JournalPubDate'; end
                                    
                                    data{end+1} = ArticleTitle;
                                    colnames{end+1} = 'ArticleTitle';
                                    
                                    if ~isempty(ArticleDate); data{end+1} = ArticleDate; colnames{end+1} = 'ArticleDate'; end
                                    
                                    if ~isempty(Pagination); data{end+1} = Pagination; colnames{end+1} = 'Pagination'; end
                                    
                                    if ~isempty(Abstract_empty); data{end+1} = Abstract_empty; colnames{end+1} = 'AbstractEmpty'; end
                                    
                                    if ~isempty(Abstract_unassigned); data{end+1} = Abstract_unassigned; colnames{end+1} = 'AbstractUnassigned'; end
                                    
                                    if ~isempty(Abstract_unlabelled); data{end+1} = Abstract_unlabelled; colnames{end+1} = 'AbstractUnlabelled'; end
                                    
                                    if ~isempty(Abstract_complete); data{end+1} = Abstract_complete; colnames{end+1} = 'AbstractComplete'; end
                                    
                                    if ~isempty(Abstract_introduction); data{end+1} = Abstract_introduction; colnames{end+1} = 'AbstractIntroduction'; end
                                    
                                    if ~isempty(Abstract_background); data{end+1} = Abstract_background; colnames{end+1} = 'AbstractBackground'; end
                                    
                                    if ~isempty(Abstract_purpose); data{end+1} = Abstract_purpose; colnames{end+1} = 'AbstractPurpose'; end
                                    
                                    if ~isempty(Abstract_objective); data{end+1} = Abstract_objective; colnames{end+1} = 'AbstractObjective'; end
                                    
                                    if ~isempty(Abstract_aim); data{end+1} = Abstract_aim; colnames{end+1} = 'AbstractAIM'; end
                                    
                                    if ~isempty(Abstract_methods); data{end+1} = Abstract_methods; colnames{end+1} = 'AbstractMethods'; end
                                    
                                    if ~isempty(Abstract_search_methods); data{end+1} = Abstract_search_methods; colnames{end+1} = 'AbstractSearchMethods'; end
                                    
                                    if ~isempty(Abstract_selection_criteria); data{end+1} = Abstract_selection_criteria; colnames{end+1} = 'AbstractSelectionCriteria'; end
                                    
                                    if ~isempty(Abstract_material); data{end+1} = Abstract_material; colnames{end+1} = 'AbstractMaterial'; end
                                    
                                    if ~isempty(Abstract_material_and_methods); data{end+1} = Abstract_material_and_methods; colnames{end+1} = 'AbstractMaterialAndMethods'; end
                                    
                                    if ~isempty(Abstract_data_collection_and_analysis); data{end+1} = Abstract_data_collection_and_analysis; colnames{end+1} = 'AbstractDataCollectionAndAnalysis'; end
                                    
                                    if ~isempty(Abstract_results); data{end+1} = Abstract_results; colnames{end+1} = 'AbstractResults'; end
                                    
                                    if ~isempty(Abstract_main_results); data{end+1} = Abstract_main_results; colnames{end+1} = 'AbstractMainResults'; end
                                    
                                    if ~isempty(Abstract_conclusion); data{end+1} = Abstract_conclusion; colnames{end+1} = 'AbstractConclusion'; end
                                    
                                    if ~isempty(Abstract_authors_conclusions); data{end+1} = Abstract_authors_conclusions; colnames{end+1} = 'AbstractAuthorsConclusions'; end
                                    
                                    if ~isempty(Abstract_importance); data{end+1} = Abstract_importance; colnames{end+1} = 'AbstractImportance'; end
                                    
                                    if ~isempty(Abstract_design_and_setting); data{end+1} = Abstract_design_and_setting; colnames{end+1} = 'AbstractDesignAndSetting'; end
                                    
                                    if ~isempty(Abstract_exposure); data{end+1} = Abstract_exposure; colnames{end+1} = 'AbstractExposure'; end
                                    
                                    if ~isempty(Abstract_main_outcomes_and_measures); data{end+1} = Abstract_main_outcomes_and_measures; colnames{end+1} = 'AbstractMainOutcomesAndMeasures'; end
                                    
                                    if ~isempty(Abstract_conclusions_and_relevance); data{end+1} = Abstract_conclusions_and_relevance; colnames{end+1} = 'AbstractConclusionsAndRelevance'; end
                                    
                                    if ~isempty(Abstract_findings); data{end+1} = Abstract_findings; colnames{end+1} = 'AbstractFindings'; end
                                    
                                    if ~isempty(Abstract_options); data{end+1} = Abstract_options; colnames{end+1} = 'AbstractOptions'; end
                                    
                                    if ~isempty(Abstract_implications); data{end+1} = Abstract_implications; colnames{end+1} = 'AbstractImplications'; end
                                    
                                    if ~isempty(Abstract_outcomes); data{end+1} = Abstract_outcomes; colnames{end+1} = 'AbstractOutcomes'; end
                                    
                                    if ~isempty(Abstract_interpretation); data{end+1} = Abstract_interpretation; colnames{end+1} = 'AbstractInterpretation'; end
                                    
                                    if ~isempty(Abstract_funding); data{end+1} = Abstract_funding; colnames{end+1} = 'AbstractFunding'; end
                                    
                                    if ~isempty(Abstract_evidence); data{end+1} = Abstract_evidence; colnames{end+1} = 'AbstractEvidence'; end
                                    
                                    if ~isempty(Abstract_context); data{end+1} = Abstract_context; colnames{end+1} = 'AbstractContext'; end
                                       
                                    if ~isempty(Abstract_values); data{end+1} = Abstract_values; colnames{end+1} = 'AbstractValues'; end
                                     
                                    if ~isempty(Abstract_design); data{end+1} = Abstract_design; colnames{end+1} = 'AbstractDesign'; end
                                     
                                    if ~isempty(Abstract_setting); data{end+1} = Abstract_setting; colnames{end+1} = 'AbstractSetting'; end
                                    
                                    if ~isempty(Abstract_recommendations); data{end+1} = Abstract_recommendations; colnames{end+1} = 'AbstractRecommendations'; end
                                    
                                    if ~isempty(Abstract_main_outcome_measures); data{end+1} = Abstract_main_outcome_measures; colnames{end+1} = 'AbstractMainOutcomeMeasures'; end
                                    
                                    if ~isempty(Abstract_statement_of_significance); data{end+1} = Abstract_statement_of_significance; colnames{end+1} = 'AbstractStatementOfSignificance'; end
                                    
                                    if ~isempty(Abstract_patients_and_methods); data{end+1} = Abstract_patients_and_methods; colnames{end+1} = 'AbstractPatientsAndMethods'; end
                                    
                                    if ~isempty(Abstract_subjects_and_methods); data{end+1} = Abstract_subjects_and_methods; colnames{end+1} = 'AbstractSubjectsAndMethods'; end
                                    
                                    if ~isempty(Abstract_background_and_methods); data{end+1} = Abstract_background_and_methods; colnames{end+1} = 'AbstractBackgroundAndMethods'; end
                                    
                                    if ~isempty(Abstract_methods_and_results); data{end+1} = Abstract_methods_and_results; colnames{end+1} = 'AbstractMethodsAndResults'; end
                                    
                                    if ~isempty(Abstract_benefits_costs_and_harms); data{end+1} = Abstract_benefits_costs_and_harms; colnames{end+1} = 'AbstractBenefitsCostsAndHarms'; end
                                    
                                    if ~isempty(Abstract_validation); data{end+1} = Abstract_validation; colnames{end+1} = 'AbstractValidation'; end
                                    
                                    if ~isempty(Abstract_sponsor); data{end+1} = Abstract_sponsor; colnames{end+1} = 'AbstractSponsor'; end
                                    
                                    data{end+1} = WHOLE_Abstract; colnames{end+1} = 'WHOLE_Abstract';
                                 
                                    if iscell(Language)
                                        
                                        x = cell.empty;
                                        for i=1:length(Language)
                                            if i==1; x = Language{1}; else x = strcat(x,'-',Language{i}); end
                                        end
                                        Language = x;
                                        
                                    else
                                        
                                        data{end+1} = Language;
                                        colnames{end+1} = 'Language';
                                        
                                    end
                                    
                                    if ~isempty(articleID_pubmed); data{end+1} = articleID_pubmed; colnames{end+1} = 'articleIDPUBMED'; end
                                    
                                    if ~isempty(articleID_doi); data{end+1} = articleID_doi; colnames{end+1} = 'articleIDDOI'; end
                                    
                                    if ~isempty(articleID_pii); data{end+1} = articleID_pii; colnames{end+1} = 'articleIDPII'; end
                                    
                                    if ~isempty(articleID_pmc); data{end+1} = articleID_pmc; colnames{end+1} = 'articleIDPMC'; end
                                    
                                    if ~isempty(articleID_mid); data{end+1} = articleID_mid; colnames{end+1} = 'articleIDMID'; end
                                    
                                    if ~isempty(articleID_pmcid); data{end+1} = articleID_pmcid; colnames{end+1} = 'articleIDPMCID'; end
                                    
                                    data{end+1} = PublicationType;
                                    colnames{end+1} = 'PublicationType';
                                                                      
                                    %%% if iPub == 1; save('tmp.mat','data','colnames'); end
                                    
                                   try

                                       %% conn = checkConn(conn);
                                       
                                        fastinsert(conn,table_name,colnames,data);
                                        clear colnames data

                                        my_commit(conn);

                                        disp(strcat('Commit:',int2str(iYear),'-',int2str(iFile)));

                                        fprintf(fileID,'%s\n',strcat('Commit:',table_name,':',int2str(iYear),'-',int2str(iFile)));


                                    catch ME

                                        rollback(conn);

                                        msg = getReport(ME);

                                        disp(msg);

                                        disp(strcat('ErrorOnFastInsert:',table_name,':',int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));

                                        disp(msg);
                                        
                                        fprintf(fileID,'%s\n',strcat('ErrorOnFastInsert:',table_name,':',int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));

                                        fprintf(fileID,'%s\n',msg);
                                        
                                        disp(strcat('Rollback:',int2str(iYear),'-',int2str(iFile)));

                                        fprintf(fileID,'%s\n',strcat('Rollback:',int2str(iYear),'-',int2str(iFile)));
                                        
                                        break 

                                   end
                                    
                                    %% conn = checkConn(conn);
                                    
                                    select_check_ArticleID = strcat('SELECT idPubmedArticle FROM pubmedarticle WHERE pubmedarticle.PMID = ',PMID);

%                                     curs = exec(conn,select_check_ArticleID);
%                                     curs = fetch(curs);

                                    curs = my_exec(conn,select_check_ArticleID);

                                    idPubmedArticle =  curs.Data{1};
                                    
                                    close(curs);
                                    
                                    if ~isempty(All_Authors)
                                       
                                        nAuthors = length(All_Authors);
                                        
                                        for iAuthor=1:nAuthors
                                           
                                            %%% CHECK IF AUTHORS EXIST
                                            
                                            LastName = All_Authors(iAuthor).LastName;
                                            ForeName = All_Authors(iAuthor).ForeName;
                                            Initials = All_Authors(iAuthor).Initials;
                                            
                                            if isempty(LastName); LastName = '(empty)'; end
                                            if isempty(ForeName); ForeName = '(empty)'; end
                                            if isempty(Initials); Initials = '(empty)'; end
                                            
                                            %% conn = checkConn(conn);
                                            
                                            select_check_Author = strcat('SELECT COUNT(*) FROM authors WHERE authors.LastName = "',LastName,'" AND authors.ForeName = "',ForeName,'"');
                                        
%                                             curs = exec(conn,select_check_Author);
%                                             curs = fetch(curs);

                                            curs = my_exec(conn,select_check_Author);

                                            if curs.Data{1} == 0; Author_exist = 0; else Author_exist = 1; end
                                            
                                            close(curs);
                                            
                                            %%% GET idAuthors
                                            
                                            if Author_exist
                                            
                                                LastName = All_Authors(iAuthor).LastName;
                                                ForeName = All_Authors(iAuthor).ForeName;
                                                Initials = All_Authors(iAuthor).Initials;
                                            
                                                if isempty(LastName); LastName = '(empty)'; end
                                                if isempty(ForeName); ForeName = '(empty)'; end
                                                if isempty(Initials); Initials = '(empty)'; end
                                            
                                                %% conn = checkConn(conn);
                                                
                                                select_check_Author = strcat('SELECT idAuthors FROM authors WHERE authors.LastName = "',LastName,'" AND authors.ForeName = "',ForeName,'"');
                                        
%                                                 curs = exec(conn,select_check_Author);
%                                                 curs = fetch(curs);

                                                curs = my_exec(conn,select_check_Author);

                                                idAuthors = curs.Data{1};
                                                
                                                close(curs);
                                                
                                            else
                                                
                                                idAuthors = cell.empty;
                                                
                                            end
                                                
                                            if isempty(idAuthors)
                                                
                                                LastName = All_Authors(iAuthor).LastName;
                                                ForeName = All_Authors(iAuthor).ForeName;
                                                Initials = All_Authors(iAuthor).Initials;
                                            
                                                if isempty(LastName); LastName = '(empty)'; end
                                                if isempty(ForeName); ForeName = '(empty)'; end
                                                if isempty(Initials); Initials = '(empty)'; end
                                            
                                                %%% INSERT AUTHOR DATA
                                                %%% FIRST
                                                table_name = 'authors';
                                                data = {LastName,ForeName,Initials};
                                                colnames = {'LastName','ForeName','Initials'};
                                                
                                               try

                                                   %% conn = checkConn(conn);

                                                    fastinsert(conn,table_name,colnames,data);
                                                    clear colnames data

                                                    my_commit(conn);

                                                    disp(strcat('Commit:',int2str(iYear),'-',int2str(iFile)));

                                                    fprintf(fileID,'%s\n',strcat('Commit:',table_name,':',int2str(iYear),'-',int2str(iFile)));


                                                catch ME

                                                    rollback(conn);

                                                    msg = getReport(ME);

                                                    disp(strcat('ErrorOnFastInsert:',table_name,':',int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));

                                                    disp(msg);
                                                    
                                                    fprintf(fileID,'%s\n',strcat('ErrorOnFastInsert:',table_name,':',int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));

                                                    fprintf(fileID,'%s\n',msg);
                                                    
                                                    disp(strcat('Rollback:',int2str(iYear),'-',int2str(iFile)));

                                                    fprintf(fileID,'%s\n',strcat('Rollback:',int2str(iYear),'-',int2str(iFile)));
                                                    
                                                    break 

                                               end
                                                
                                                %% conn = checkConn(conn);
                                                
                                                select_check_Author = strcat('SELECT idAuthors FROM authors WHERE authors.LastName = "',LastName,'" AND authors.ForeName = "',ForeName,'"');
                                        
%                                                 curs = exec(conn,select_check_Author);
%                                                 curs = fetch(curs);

                                                curs = my_exec(conn,select_check_Author);

                                                idAuthors = curs.Data{1};
                                                
                                                close(curs);
                                                
                                            end
                                                
                                            %%% INSERT ARTICLE DATA
                                            
                                            table_name = 'authorsperpublication';
                                            
                                            idAuthors_cell = idAuthors;
                                            idPubmedArticle_cell = idPubmedArticle;
                                            
                                            fields_All_Authors = fieldnames(All_Authors(iAuthor));
                                            
                                            AffiliationInfo_cell = cell.empty;
                                            
                                            if sum(strcmp('AffiliationInfo',fields_All_Authors)) ~= 0
                                                
                                                if isstruct(All_Authors(iAuthor).AffiliationInfo)
                                                
                                                    fields_AffiliationInfo = fieldnames(All_Authors(iAuthor).AffiliationInfo);
                                                
                                                    if sum(strcmp('Affiliation',fields_AffiliationInfo)) ~= 0
                                                
                                                        AffiliationInfo_cell{1} = All_Authors(iAuthor).AffiliationInfo.Affiliation;
                                                    
                                                    end
                                                    
                                                end

                                            end
                                            
                                            data{1} = idAuthors_cell;
                                            colnames{1} = 'AuthorsID';
                                            
                                            data{end+1} = idPubmedArticle_cell;
                                            colnames{end+1} = 'PubmedArticleID';
                                            
                                            if ~isempty(AffiliationInfo_cell); data{end+1} = AffiliationInfo_cell{1}; colnames{end+1} = 'Affiliation'; end
                                            
                                             try

                                                 %% conn = checkConn(conn);
                                                 
                                                    fastinsert(conn,table_name,colnames,data);
                                                    clear colnames data

                                                    my_commit(conn);

                                                    disp(strcat('Commit:',int2str(iYear),'-',int2str(iFile)));

                                                    fprintf(fileID,'%s\n',strcat('Commit:',table_name,':',int2str(iYear),'-',int2str(iFile)));


                                             catch ME

                                                    rollback(conn);

                                                    msg = getReport(ME);

                                                    disp(strcat('ErrorOnFastInsert:',table_name,':',int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));

                                                    disp(msg);
                                                    
                                                    fprintf(fileID,'%s\n',strcat('ErrorOnFastInsert:',table_name,':',int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));

                                                    fprintf(fileID,'%s\n',msg);
                                                    
                                                    disp(strcat('Rollback:',int2str(iYear),'-',int2str(iFile)));

                                                    fprintf(fileID,'%s\n',strcat('Rollback:',int2str(iYear),'-',int2str(iFile)));
                                                    
                                                    break 

                                            end
                                            
                                        end
                                        
                                    end
                                    
                                    if ~isempty(All_Keywords)
                                       
                                        nKeywords = length(All_Keywords);
                                        
                                        for iKeyword=1:nKeywords
                                           
                                            %%% CHECK IF KEYWORD EXIST
                                            
                                            %% conn = checkConn(conn);
                                            
                                            select_check_Keyword = strcat('SELECT COUNT(*) FROM keywords WHERE keywords.Keyword = "',All_Keywords(iKeyword).CONTENT,'"');
                                        
%                                             curs = exec(conn,select_check_Keyword);
%                                             curs = fetch(curs);
                                            

                                            curs = my_exec(conn,select_check_Keyword);
                                            
                                            if curs.Data{1} == 0; Keyword_exist = 0; else Keyword_exist = 1; end
                                            
                                            close(curs);
                                            
                                             %%% GET idAuthors
                                            
                                            if Keyword_exist
                                                
                                                %% conn = checkConn(conn);
                                                
                                                select_check_Keyword = strcat('SELECT idKeywords FROM keywords WHERE keywords.Keyword = "',All_Keywords(iKeyword).CONTENT,'"');
                                        
                                                % curs = exec(conn,select_check_Keyword);
                                                % curs = fetch(curs);
                                                
                                                curs = my_exec(conn,select_check_Keyword);
                                                
                                                
                                                idKeywords = curs.Data{1};
                                                
                                                close(curs);
                                                
                                            else
                                                
                                                idKeywords = cell.empty;
                                                
                                            end
                                            
                                             if isempty(idKeywords)
                                                
                                                %%% INSERT KEYWORD DATA
                                                %%% FIRST
                                                table_name = 'keywords';
                                                data = {All_Keywords(iKeyword).CONTENT};
                                                colnames = {'Keyword'};
                                    
                                                try

                                                    %% conn = checkConn(conn);
                                                    
                                                    fastinsert(conn,table_name,colnames,data);
                                                    clear colnames data

                                                    my_commit(conn);

                                                    disp(strcat('Commit:',int2str(iYear),'-',int2str(iFile)));

                                                    fprintf(fileID,'%s\n',strcat('Commit:',table_name,':',int2str(iYear),'-',int2str(iFile)));


                                               catch ME

                                                    rollback(conn);

                                                    msg = getReport(ME);

                                                    disp(strcat('ErrorOnFastInsert:',table_name,':',int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));

                                                    disp(msg);
                                                    
                                                    fprintf(fileID,'%s\n',strcat('ErrorOnFastInsert:',table_name,':',int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));

                                                    fprintf(fileID,'%s\n',msg);
                                                    
                                                    disp(strcat('Rollback:',int2str(iYear),'-',int2str(iFile)));

                                                    fprintf(fileID,'%s\n',strcat('Rollback:',int2str(iYear),'-',int2str(iFile)));
                                                    
                                                    break

                                                end
                                                
                                                select_check_Keyword = strcat('SELECT idKeywords FROM keywords WHERE keywords.Keyword = "',All_Keywords(iKeyword).CONTENT,'"');
                                        
                                                %curs = exec(conn,select_check_Keyword);
                                                %curs = fetch(curs);
                                                
                                                curs = my_exec(conn,select_check_Keyword);
                                                
                                                
                                                idKeywords = curs.Data{1};
                                                
                                                close(curs);
                                                
                                             end
                                             
                                            %%% INSERT ARTICLE DATA
                                            table_name = 'keywordsperpublication';
                                            data = {idKeywords,idPubmedArticle};
                                            colnames = {'KeywordsID','PubmedArticleID'};

                                            try

                                                %% conn = checkConn(conn);
                                                
                                                fastinsert(conn,table_name,colnames,data);
                                                clear colnames data

                                                my_commit(conn);

                                                disp(strcat('Commit:',int2str(iYear),'-',int2str(iFile)));

                                                fprintf(fileID,'%s\n',strcat('Commit:',table_name,':',int2str(iYear),'-',int2str(iFile)));


                                            catch ME

                                                rollback(conn);

                                                msg = getReport(ME);

                                                disp(strcat('ErrorOnFastInsert:',table_name,':',int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));

                                                disp(msg);
                                                
                                                fprintf(fileID,'%s\n',strcat('ErrorOnFastInsert:',table_name,':',int2str(iYear),'-',int2str(iFile),'-',int2str(iPub)));

                                                fprintf(fileID,'%s\n',msg);
                                                
                                                disp(strcat('Rollback:',int2str(iYear),'-',int2str(iFile)));

                                                fprintf(fileID,'%s\n',strcat('Rollback:',int2str(iYear),'-',int2str(iFile)));
                                                
                                                break

                                           end
                                            
                                        end
                                        
                                    end
                                    
                                end
                              
                          end
                          
                      end
 
                      close(conn);
                      
                      fclose(fileID);
           
                  end
                  
              end
              
           end
          
       end
       
       end
        
    end

% end

end

end

function conn = checkConn(conn)

isONLINE = false;

wait_database_time = 120;

firstTime = 1;

while ~(isONLINE)
    
    if ~(firstTime)
        
        pause(wait_database_time); 
        
        disp('...waiting for database');
    
        conn = database('references_papers','joaodornas','Senh@1111','Vendor','MySQL','Server','141.44.42.20','PortNumber',7979);
    
    end
    
    try 

        info = ping(conn);
        
        disp('...ping to database done');

        expected = 'joaodornas@INDIREA-SERVER';

        currentUser = info.CurrentUserName;

        if strfind(currentUser,expected)

            isONLINE = true;

        else 

            isONLINE = false;
            
            clear conn

        end

    catch ME
        
        msg = getReport(ME);
        
        disp(msg);
        
        isONLINE = false;
        
        clear conn

    end
    
    firstTime = 0;

end


end

function my_commit(conn)

isONLINE = false;

wait_database_time = 120;

firstTime = 1;

while ~(isONLINE)
    
    if ~(firstTime)
        
        pause(wait_database_time); 
        
        disp('...waiting for database');

    end
    
    try 

        commit(conn);
        
        disp('...commit done');

        isONLINE = true;

    catch ME
        
        msg = getReport(ME);
        
        disp(msg);
        
        isONLINE = false;
        
    end
    
    firstTime = 0;

end

end

function curs = my_exec(conn,selectStatement)

isONLINE = false;

wait_database_time = 3*60;

firstTime = 1;

while ~(isONLINE)
    
    if ~(firstTime)
        
        pause(wait_database_time); 
        
        disp('...waiting for database');

    end
    
    try 

        curs = exec(conn,selectStatement);
        
        if isempty(curs.Message)
        
            curs = fetch(curs);
        
            disp('...select done');
        
            isONLINE = true;
            
        else
            
            disp(curs.Message);
            
            clear conn
            conn = database('references_papers','joaodornas','Senh@1111','Vendor','MySQL','Server','141.44.42.20','PortNumber',7979);
            
        end

    catch ME
        
        msg = getReport(ME);
        
        disp(msg);
        
        isONLINE = false;
        
    end
    
    firstTime = 0;

end

end

