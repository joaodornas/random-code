

%%% READ REVIEWS INFORMATION AND GET UNIQUE KEYWORDS

folder = '/Volumes/DORNAS/__REVIEWS/_Retina-Reviews';

files = dir(strcat(folder,'/','*.csv'));

allKeys = [];

nFiles = length(files);

for iFile=1:nFiles
   
   name = files(iFile).name;
   
   if ~strcmp(name(1:2),'._')
       
       disp(name);
    
       fid = fopen(strcat(folder,'/',name)); 

       data = textscan(fid,'%s');

       fclose(fid);

       tmp_keys = data{1};

       words_cleaned = cleanWords(tmp_keys);

       allKeys = [allKeys,words_cleaned];
       
   end
   
end

uKeys = unique(allKeys);
uKeys = uKeys';