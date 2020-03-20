
list = dir('*-json.html');

nFiles = length(list);

all_fields = cell.empty;

for iFile=1:nFiles
    
   disp(int2str(iFile));
   
   try
    
        jsonfile = loadjson(list(iFile).name);
   
        getFields = fields(jsonfile);
   
        all_fields = [all_fields; getFields];
        
   catch ME
       
       msg = getReport(ME);
       
       disp(msg);
       
   end
    
end

u_fields = unique(all_fields);

disp('FINISHED !!');

save('all_fields.mat','all_fields','u_fields');