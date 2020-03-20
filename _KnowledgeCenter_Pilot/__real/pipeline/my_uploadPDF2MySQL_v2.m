
% javaaddpath('/Users/joaodornas/Dropbox (joaodornas)/_Research/_CODES/reference-database/Reference-Database/DB_drivers/mysql-connector-java-5.1.40-bin.jar');

project_label = 'Retina';
% source_label = 'DOWNLOADED'; %%% 1. DOWNLOADED 2. PAPERSFORMAC
source_label = 'PAPERSFORMAC'; %%% 1. DOWNLOADED 2. PAPERSFORMAC

if strcmp(project_label,'Retina')
    
    query_id = 1;
    
    mainFilterLabel = 'HEALTHY';
    
    main_filter_version = 1;
    
    secondFilter = 'PRIMATES';
    
    second_filter_version = 4;
    
    database_label = 'retina';
    
    conn = database(strcat('docs_papers_',database_label),'root','senha1111','Vendor','MySQL','Server','127.0.0.1','PortNumber',3306);
    
    main_doc_folder = strcat('/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/__data/DOCS_PAPERS','/'); 
    
    main_project_folder = strcat('/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/Visual-Attention-Awareness/2. Retina Model/_ArticleTitle_v',int2str(second_filter_version),'/');
    
    %%% queryName = " RETINA:1:HEALTHY:1:PRIMATES:4 "
    
    queryName = strcat(upper(project_label),':',int2str(query_id),':',mainFilterLabel,':',int2str(main_filter_version),':',secondFilter,':',int2str(second_filter_version));   
    
    HTML_SOURCE{1} = 'ADOBE';
    HTML_SOURCE{2} = 'PDFBOX';
    
end

selectFieldID = sprintf('SELECT idqueries FROM queries WHERE queryName = ''%s''',queryName);

curs = exec(conn,selectFieldID);
curs = fetch(curs);

if curs.Data{1} == 0; queryID_exist = 0; else queryID_exist = 1; end

close(curs);

if queryID_exist

    queryID = curs.Data{1};

    pdf_source_folder = strcat(main_doc_folder,source_label,'/','PDF');

    load(strcat(main_project_folder,'__selectedPMIDs.mat'));

    PMIDs = [];

    eval(sprintf('PMIDs = %s;',source_label));

    nPMID = length(PMIDs);

    for iPMID=1:nPMID

        try 

            PMID = PMIDs(iPMID);

        %%%% LOAD DOCS

            clear data colnames

            table_name = 'docs';

            selectFieldID = strcat('SELECT iddoc FROM docs WHERE docs.PMID = ',int2str(PMID),' AND SOURCE = "',source_label,'"');

            curs = exec(conn,selectFieldID);
            curs = fetch(curs);
            
            response = curs.Data(1);

            if sum(response{1} == 0) || strcmp(response{1},'No Data'); docID_exist = 0; else docID_exist = 1; end

            close(curs);

            if docID_exist

                disp(strcat('PMID exists:',int2str(PMID)));

                iddoc = curs.Data(1);

                %%% DO DOC PER QUERY

                clear data colnames

                table_name = 'docsperqueries';

                data{1} = iddoc{1};
                colnames{1} = 'iddoc';

                data{2} = queryID;
                colnames{2} = 'idquery';

                fastinsert(conn,table_name,colnames,data);

            else

                disp(strcat('PMID will be uploaded:',int2str(PMID)));

                %%% ADD DOC

                data{1} = PMID;
                colnames{1} = 'PMID';

                data{2} = source_label;
                colnames{2} = 'SOURCE';

                fastinsert(conn,table_name,colnames,data);

                selectFieldID = strcat('SELECT iddoc FROM docs WHERE docs.PMID = ',int2str(PMID),' AND SOURCE = "',source_label,'"');

                curs = exec(conn,selectFieldID);
                curs = fetch(curs);

                iddoc = curs.Data{1};

                %%% DO DOC PER QUERY

                clear data colnames

                table_name = 'docsperqueries';

                data{1} = iddoc;
                colnames{1} = 'iddoc';

                data{2} = queryID;
                colnames{2} = 'idquery';

                fastinsert(conn,table_name,colnames,data);

                %%%% LOAD METADATA

                meta_source_folder = strcat(main_doc_folder,source_label,'/','METADATA-TIKA');

                clear data colnames

                table_name = 'metadataperdoc';

                metadata = loadjson(strcat(meta_source_folder,'/',int2str(PMID),'-metadata-json.html'));

                theseFields = fields(metadata);

                for iField=1:length(theseFields)

                   thisField = theseFields(iField);
                   thisField = thisField{1};

                   value = [];

                   eval(strcat('value=metadata.',thisField,';'));

                   selectFieldID = strcat('SELECT MetaDataFieldsID FROM metadatafields WHERE MetaDataFieldsName = "',thisField,'"');

                   curs = exec(conn,selectFieldID);
                   curs = fetch(curs);

                   if curs.Data{1} == 0; MetaID_exist = 0; else MetaID_exist = 1; end

                   close(curs);

                   if MetaID_exist

                       clear data colnames

                       data{1} = iddoc;
                       colnames{1} = 'iddoc';

                       data{2} = curs.Data{1};
                       colnames{2} = 'idmetadata';

                       if iscell(value); new_value = value(1); else new_value = value; end

                       data{3} = new_value;
                       colnames{3} = 'metadatavalue';

                       fastinsert(conn,table_name,colnames,data);

                   else

                       break

                   end

                end

                %%%% LOAD HTML - SOURCE 1:ADOBE

                htmltable_name = 'htmltextperdocs';
                htmltable_name2 = 'htmlimagesperdocs';

                html_source_folder = strcat(main_doc_folder,source_label,'/','HTML-',HTML_SOURCE{1});

                out_files = dir(strcat(html_source_folder,'/',int2str(PMID),'.html'));

                in_files = dir(strcat(html_source_folder,'/',int2str(PMID),'/','*.htm'));

                whole_html = [];
                for iFile=1:length(out_files)

                    fid = fopen(strcat(html_source_folder,'/',out_files(iFile).name),'r','n','utf-8');
                    f = textscan(fid,'%s ');
                    ff = f{1};
                    whole_html = [whole_html;ff];

                    fclose(fid);

                end

                for iFile=1:length(in_files)

                    name = strcat('part',int2str(iFile),'.htm');

                    fid = fopen(strcat(html_source_folder,'/',int2str(PMID),'/',name),'r','n','utf-8');
                    f = textscan(fid,'%s ');
                    ff = f{1};
                    whole_html = [whole_html;ff];

                    fclose(fid);

                end

                if ~isempty(whole_html); whole_html = strjoin(whole_html,' '); end

                clear data colnames

                data{1} = iddoc;
                colnames{1} = 'iddoc';

                data{2} = whole_html;
                colnames{2} = 'html';

                data{3} = HTML_SOURCE{1};
                colnames{3} = 'htmlsource';

                fastinsert(conn,htmltable_name,colnames,data);

                %%%% LOAD HTML - SOURCE 2:PDFBOX

                htmltable_name = 'htmltextperdocs';

                html_source_folder = strcat(main_doc_folder,source_label,'/','HTML-',HTML_SOURCE{2});

                files = dir(strcat(html_source_folder,'/',int2str(PMID),'*.txt'));

                for iFile=1:length(files)

                    clear data colnames

                    whole_html = [];

                    fid = fopen(strcat(html_source_folder,'/',files(iFile).name),'r','n','utf-8');
                    f = textscan(fid,'%s ');
                    whole_html = f{1};

                    fclose(fid);

                    whole_html = strjoin(whole_html,' ');

                    data{1} = iddoc;
                    colnames{1} = 'iddoc';

                    data{2} = whole_html;
                    colnames{2} = 'html';

                    data{3} = HTML_SOURCE{2};
                    colnames{3} = 'htmlsource';

                    fastinsert(conn,htmltable_name,colnames,data);

                end


            end

        catch ME

            msg = getReport(ME);

            disp(strcat('ERROR ON PMID:',int2str(PMID)));

            disp(msg);

            fid = fopen(strcat(int2str(PMID),'-error.txt'),'w');
            fclose(fid);

        end

    end
    
else

    disp('Query does not exist !!');

end


