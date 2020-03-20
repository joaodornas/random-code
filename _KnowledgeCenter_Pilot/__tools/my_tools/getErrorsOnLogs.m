
files = dir('*.txt');

nFiles = length(files);

fid_error = fopen('error-retina-1.txt','a+');

for iFile=1:nFiles
    
   fid = fopen(files(iFile).name,'r');
   
   save_line = 0;
   stop_save = 0;
   
   rline = fgets(fid);
   while ischar(rline)
  
       if strfind(lower(rline),'error')

          save_line = 1;
           
       elseif strfind(lower(rline),'rollback')
           
           stop_save = 1;
           
       end
       
       if save_line
           
           fprintf(fid_error,'%s\n',rline);
           
           if stop_save == 1;
               
               save_line = 0;
               stop_save = 0;
               
           end
           
       end
       
       rline = fgets(fid);
       
   end
   
   fclose(fid);
   
end

fclose(fid_error);
   