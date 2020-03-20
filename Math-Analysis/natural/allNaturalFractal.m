function allNaturalFractal(start_registro,end_registro,OS)

registro = importdata('memoryBackwardProtocols.txt');

if strcmp(OS,'Mac')
    
    folderPath = strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fractal/');
    bar = '/';
    
elseif strcmp(OS,'Win')
    
    folderPath = strcat('Z:\DATA\Forward-Backward\Power-Law\Fractal\');
    bar = '\';
    
end

for r=start_registro:end_registro
    
   disp(registro{r}); 
    
   protocol = registro{r}(2:end-4);
   
   mainPathDong = strcat(folderPath,'Dong',bar);
   mainPathDornas = strcat(folderPath,'Dornas',bar);
   
   analysisFractalDongDornas(protocol,bar,'dong',mainPathDong);
   
   analysisFractalDongDornas(protocol,bar,'dornas',mainPathDornas);
   
   
end


end

