
tFile_int = [1000 2000 3000 4000 5000 6000 7000];

for iMen=tFile_int
    
    load(strcat('Mendeley_',int2str(iMen),'.mat')); %%% Mendeley

    nPDF = length(file);

    count_without_pmid = 0;

    % wrong_ID = 0;
    for iFile=1:nPDF

       this_file = file(iFile); 
       this_pmid = pmid(iFile);

       %%% PAPERS FOR MAC

    %    if length(this_pmid{:}) > 1 && length(this_file{:}) > 1
    %       
    %        idx_start = strfind(this_file,'/Users');
    %        idx_end = strfind(this_file,'.pdf:application');
    %        vec = this_file{:};
    %        
    %        whole_file_path = vec(idx_start{:}:idx_end{:}+3);
    %        
    %        idx_pmid_start = strfind(this_pmid,'{');
    %        idx_pmid_end = strfind(this_pmid,'}');
    %        vec = this_pmid{:};
    % 
    %        PMID = str2num(vec(idx_pmid_start{:}+1:idx_pmid_end{:}-1));
    %        
    %        system(sprintf('cp "%s" pdfs/%s.pdf',whole_file_path,int2str(PMID)));
    %        
    %    end

    %%% READ CUBE
    % 
    %    if length(this_pmid{:}) > 1 && length(this_file{:}) > 1
    %        
    %        whole_file_path = this_file{:};
    % 
    %        PMID = str2num(this_pmid{:});
    %        
    %        system(sprintf('mv "%s" pdfs/%s.pdf',whole_file_path,int2str(PMID)));
    %        
    %    end

       %%% MENDELEY

       if length(this_pmid{:}) > 1 && length(this_file{:}) > 1

           whole_file_path = this_file{:};

           PMID = str2num(this_pmid{:});

           system(sprintf('mv "/%s" Mendeley/%s.pdf',whole_file_path,int2str(PMID)));

       end

    %%% PAPERS FOR MAC
    %    if length(this_pmid{:}) <= 1 && length(this_file{:}) > 1
    %        
    %        count_without_pmid = count_without_pmid + 1;
    %        
    %    end

    %    if length(this_pmid{:}) <= 1 && length(this_file{:}) > 1
    %        
    %        wrong_ID = wrong_ID + 1;
    %       
    %        idx_start = strfind(this_file,'/Users');
    %        idx_end = strfind(this_file,'.pdf:application');
    %        vec = this_file{:};
    %        
    %        whole_file_path = vec(idx_start{:}:idx_end{:}+3);
    %        
    %        system(sprintf('cp "%s" pdfs/%s.pdf',whole_file_path,int2str(wrong_ID)));
    %        
    %    end

    end

end
