
% tFile = 'all_papers_from_PAPERS_for_Mac.bib';
% tFile = 'readcube_export.bib'; %%% ReadCube

%%% MENDELEY
tFile_int = [1000 2000 3000 4000 5000 6000 7000];

for iMen=tFile_int
    
    tFile = strcat('Mendeley_',int2str(iMen),'.bib'); %%% Mendeley

    rBib = bibload(tFile); 

    nBibs = length(bib.cBib);

    for iBib=1:nBibs

        whole_string = bib.cBib{iBib};

        %%% PAPERS FOR MAC
    %     idx_doi = strfind(whole_string,'doi =');
    %     idx_pmid = strfind(whole_string,'pmid = ');
    %     idx_file = strfind(whole_string,'file =');
    %     idx_url = strfind(whole_string,'url = {http');

        %%% READ CUBE
    %     idx_doi = strfind(whole_string,'doi=');
    %     idx_pmid = strfind(whole_string,'pmid=');
    %     idx_file = strfind(whole_string,'file=');
    %     idx_url = strfind(whole_string,'url={http');

        %%% MENDELEY
        idx_doi = strfind(whole_string,'doi = ');
        idx_pmid = strfind(whole_string,'pmid = ');
        idx_file = strfind(whole_string,'file = ');
        idx_url = strfind(whole_string,'url = {http');
        idx_title = strfind(whole_string,'title = {{');
        idx_booktitle = strfind(whole_string,'booktitle = {');
        
        if ~isempty(idx_booktitle) && ~isempty(idx_title); idx_title(idx_title == idx_booktitle(1) + 4) = []; end

        if ~isempty(idx_doi) && idx_doi > 0

            idx_end = strfind(whole_string,'},');

            idx_end = idx_end - idx_doi;
            idx_end(idx_end<0) = [];
            min_idx_end = min(idx_end);

            doi{iBib} = whole_string(idx_doi:(idx_doi+min_idx_end));

        else

            doi{iBib} = 0;

        end

        if ~isempty(idx_pmid) && idx_pmid > 0

            idx_end = strfind(whole_string,'},');

            idx_end = idx_end - idx_pmid;
            idx_end(idx_end<0) = [];
            min_idx_end = min(idx_end);

            %%% PAPERS FOR MAC
            % pmid{iBib} = whole_string(idx_pmid:(idx_pmid+min_idx_end));

            %%% READ CUBE
            % pmid{iBib} = whole_string((idx_pmid+6):(idx_pmid+min_idx_end-1));

            %%% MENDELEY
            pmid{iBib} = whole_string((idx_pmid+8):(idx_pmid+min_idx_end-1));

        else

            pmid{iBib} = 0;

        end

        if ~isempty(idx_file) && idx_file > 0

            idx_end = strfind(whole_string,'},');

            idx_end = idx_end - idx_file;
            idx_end(idx_end<0) = [];
            min_idx_end = min(idx_end);

            %%% PAPERS FOR MAC
            % file{iBib} = whole_string(idx_file:(idx_file+min_idx_end));

            %%% READ CUBE
            % file{iBib} = whole_string((idx_file+6):(idx_file+min_idx_end-1));

            %%% MENDELEY
            file_string = whole_string((idx_file+9):(idx_file+min_idx_end-5));
            file_string = strrep(file_string,'{','');
            file_string = strrep(file_string,'}','');
            file_string = strrep(file_string,'\','');
            file{iBib} = file_string;

        else

            file{iBib} = 0;

        end
        
        if ~isempty(idx_title) && idx_title > 0

            idx_end = strfind(whole_string,'}},');

            idx_end = idx_end - idx_title;
            idx_end(idx_end<0) = [];
            min_idx_end = min(idx_end);

            %%% MENDELEY
            title_string = whole_string((idx_title+10):(idx_title+min_idx_end));
            title_string = strrep(title_string,'{','');
            title_string = strrep(title_string,'}','');
            title_string = strrep(title_string,'\','');
            title{iBib} = title_string;

        else

            title{iBib} = 0;

        end

        if ~isempty(idx_url) && idx_url > 0

            idx_end = strfind(whole_string,'},');

            idx_end = idx_end - idx_url;
            idx_end(idx_end<0) = [];
            min_idx_end = min(idx_end);

            %%% READ CUBE
            % url{iBib} = whole_string(idx_url:(idx_url+min_idx_end));

            %%% MENDELEY
            url{iBib} = whole_string(idx_url+7:(idx_url+min_idx_end-1));

        else

            url{iBib} = 0;

        end

    end

    save(strcat('Mendeley_',int2str(iMen),'.mat'));

end