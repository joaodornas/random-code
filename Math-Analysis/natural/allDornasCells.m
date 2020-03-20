function allDornasCells(start_registro,end_registro,OS)

registro = importdata('memoryBackwardProtocols.txt');

nFrames = 300;

save = 1;

if strcmp(OS,'Mac')
    
    mainPath = strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Dornas/');
    
elseif strcmp(OS,'Win')
    
    mainPath = strcat('Z:\DATA\Forward-Backward\Power-Law\Dornas\');
    
end

%for r=1:length(registro)
for r=start_registro:end_registro
    
   disp(registro{r}); 
    
   [X, Y, R] = receptiveField(char(registro{r}));
   
   if r < 27
       
       name = registro{r};
       
       name = fliplr(name);
       
       video = str2num(name(5));
       
   else
       
       video = 3;
       
   end
   
   width = X;
   
   height = Y;
   
   radius = R;
   
   [ForwardMovieCRF, BackwardMovieCRF] = ForBackCRFmovie(char(registro{r}(2:end-4)),video,nFrames,width,height,radius,save,mainPath);
   
   dornasCorrelFFT(char(registro{r}(2:end-4)),ForwardMovieCRF,'Forward',mainPath);
   
   dornasCorrelFFT(char(registro{r}(2:end-4)),BackwardMovieCRF,'Backward',mainPath);
   
   clear ForwardMovieCRF;
   clear BackwardMovieCRF;
    
end


end

