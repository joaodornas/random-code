
%%% getMetaDataFields

all_files = dir('*.html');

nFiles = length(all_files);

all_fields = cell.empty;
all_data = cell.empty;

for iFile=1:nFiles
    
    disp(int2str(iFile));
    
    fid = fopen(all_files(iFile).name,'r');

    no_meta = 0;
    rline = fgets(fid);
    while ischar(rline)
        
        idx_dots = strfind(rline,':');

        idx_created = strfind(lower(rline),'created:');
        idx_date = strfind(lower(rline),'date:');
        idx_modified = strfind(lower(rline),'modified:');

        if ~isempty(idx_date) || ~isempty(idx_created) || ~isempty(idx_modified)

            idx = idx_dots(end-2);

        else

            if ~isempty(idx_dots)
                
                idx = idx_dots(end);
                
            else
                
                no_meta = 1;
                
            end

        end
        
        if ~no_meta
            
            field_name = rline(1:idx-1);
            field_data = rline(idx+1:end);
            
            all_fields{end+1} = field_name;
            all_data{end+1} = field_data;
            
        end

        rline = fgets(fid);
        
        no_meta = 0;

    end

    fclose(fid);

end
