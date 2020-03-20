

%%% POPULATE DOCS TABLE INFORMATION

conn = database('docs_papers_retina','root','senha1111','Vendor','MySQL','Server','127.0.0.1','PortNumber',3306);

tablename = 'docs';
colnames = {'PMID','SOURCE'};

folder = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/Visual-Attention-Awareness/2. Retina Model/_ArticleTitle_v4/__PDF-1-PapersForMac/';
SOURCE = {'PapersForMac'};
files = dir(strcat(folder,'*.pdf'));

for iPMID=1:length(files)
    
    PMID(iPMID) = str2num(files(iPMID).name(1:end-4));
    
end

for iPMID=1:length(PMID)
    
    data{1} = PMID(iPMID);
    data{2} = SOURCE;

    fastinsert(conn,tablename,colnames,data);
    
end

clear PMID

folder = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/Visual-Attention-Awareness/2. Retina Model/_ArticleTitle_v4/__PDF-2-DOWNLOADED/';
SOURCE = {'DOWNLOADED'};
files = dir(strcat(folder,'*.pdf'));

for iPMID=1:length(files)
    
    PMID(iPMID) = str2num(files(iPMID).name(1:end-4));
    
end

for iPMID=1:length(PMID)
    
    data{1} = PMID(iPMID);
    data{2} = SOURCE;

    fastinsert(conn,tablename,colnames,data);
    
end




