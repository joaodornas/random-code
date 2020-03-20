
meta_folder = '/Users/joaodornas/Downloads/__Papers-for-Mac/pdfs-with-PMID-from-MISSING/metadata/';

meta_files = dir(strcat(meta_folder,'*.html'));

nFiles = length(meta_files);

doi{1} = 0;
pii{1} = 0;
pmid{1} = 0;

doi_file{1} = 0;
pii_file{1} = 0;
pmid_file{1} = 0;

for iFile=1:nFiles
    
    disp(int2str(iFile));
   
    fid = fopen(strcat(meta_folder,meta_files(iFile).name),'r');
    whole = fscanf(fid,'%s');
    fclose(fid);
    
    idx_doi = strfind(whole,'doi');
    idx_pii = strfind(whole,'PII');
    idx_pmid = strfind(whole,'PMID');
    idx_points = strfind(whole,',');
   
    if ~isempty(idx_doi)
        
        idx_point_doi = idx_points - idx_doi(1);
        idx_point_doi(idx_point_doi<0) = [];
        idx_point_doi = min(idx_point_doi);
        
        doi{end+1} = whole((idx_doi(1)+4):(idx_doi(1)+idx_point_doi-2));
        
        doi_file{end+1} = meta_files(iFile).name;
        
    end
    
    if ~isempty(idx_pii)
        
        idx_point_pii = idx_points - idx_pii(1);
        idx_point_pii(idx_point_pii<0) = [];
        idx_point_pii = min(idx_point_pii);
        
        pii{end+1} = whole((idx_pii(1)+4):(idx_pii(1)+idx_point_pii-2));
        
        pii_file{end+1} = meta_files(iFile).name;
        
    end
    
        if ~isempty(idx_pmid)
        
        idx_point_pmid = idx_points - idx_pmid(1);
        idx_point_pmid(idx_point_pmid<0) = [];
        idx_point_pmid = min(idx_point_pmid);
        
        pmid{end+1} = whole((idx_pmid(1)+4):(idx_pmid(1)+idx_point_pmid-2));
        
        pmid_file{end+1} = meta_files(iFile).name;
        
    end
    
end