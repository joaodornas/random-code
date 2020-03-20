
nMen = [1000 2000 3000 4000 5000 6000 7000];

iiFile = 0;
for iMen=nMen
    
    load(strcat('Mendeley_',int2str(iMen),'.mat'));

    for iFile=1:length(file)
        
        iiFile = iiFile + 1;

        name = file{iFile};

        if ~isempty(name)

            idx_slash = strfind(name,'/');
            
            if  ~isempty(idx_slash)
                
                [name_int,s] = str2num(name(idx_slash(end)+1:end-4));
                
            else
                
                s = 0;
                
            end

        else

            name_int = 0;

        end

        if s == 0; name_int = 0; end

        all_files(iiFile,1) = name_int;
        all_titles{iiFile,1} = title{iFile};

    end

clear title file

end