
%%% POPULATE METADATA TABLE INFORMATION

% javaaddpath('/Users/joaodornas/Dropbox (joaodornas)/_Research/_CODES/reference-database/Reference-Database/__tools/DB_drivers/mysql-connector-java-5.1.40-bin.jar');

conn = database('docs_papers_retina','root','senha1111','Vendor','MySQL','Server','127.0.0.1','PortNumber',3306);

tablename = 'metadatafields';
colnames = {'MetaDataFieldsName'};

load('uniqueMetadata.mat');

nFields = length(u);

for iField=1:nFields
    
    data{1} = u{iField};

    fastinsert(conn,tablename,colnames,data);
    
end

close(conn);
