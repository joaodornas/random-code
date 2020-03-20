function allNaturalCells(start_registro,end_registro,OS)

registro = importdata('memoryBackwardProtocols.txt');

nFrames = 300;

save = 0;

if strcmp(OS,'Mac')
    
    folderPath = strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Cena Natural/');
    bar = '/';
    
elseif strcmp(OS,'Win')
    
    folderPath = strcat('Z:\DATA\Forward-Backward\Power-Law\Cena Natural\');
    bar = '\';
    
end

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
   
   mainPath = strcat(folderPath,'Video',bar);
   
   [ForwardMovieCRF, BackwardMovieCRF] = ForBackCRFmovie(char(registro{r}(2:end-4)),video,nFrames,width,height,radius,save,mainPath);
   
    for i=1:nFrames
    
        FOR_frame(:,:) = ForwardMovieCRF(:,:,i);
        
        BACK_frame(:,:) = BackwardMovieCRF(:,:,i);
   
        reshapeframe_FOR = reshape(FOR_frame.',[],1).';
        
        reshapeframe_BACK = reshape(BACK_frame.',[],1).';
    
        FOR_reshapeframes(:,i) = reshapeframe_FOR(:); 
        
        BACK_reshapeframes(:,i) = reshapeframe_BACK(:); 
    
    end
   
   %mainPath = strcat(folderPath,'Dornas',bar);
   mainPath = strcat(folderPath,'Dornas-WO-time',bar);
   
   disp('dornasfourier-Forward');
   
   %dornasCorrelFFT(char(registro{r}(2:end-4)),ForwardMovieCRF,'Forward',mainPath);
   dornasCorrelFFTWOtime(char(registro{r}(2:end-4)),ForwardMovieCRF,'Forward',mainPath,0);
   dornasCorrelFFTWOtime(char(registro{r}(2:end-4)),ForwardMovieCRF,'Forward',mainPath,1);
   
   disp('dornasfourier-Backward');
   
   %dornasCorrelFFT(char(registro{r}(2:end-4)),BackwardMovieCRF,'Backward',mainPath);
   dornasCorrelFFTWOtime(char(registro{r}(2:end-4)),BackwardMovieCRF,'Backward',mainPath,0);
   dornasCorrelFFTWOtime(char(registro{r}(2:end-4)),BackwardMovieCRF,'Backward',mainPath,1);
   
   %mainPath = strcat(folderPath,'Dong',bar);
   mainPath = strcat(folderPath,'Dong-WO-time',bar);
   
   disp('dongfourier-Forward');
   
   %dongCorrelFFT(char(registro{r}(2:end-4)),FOR_reshapeframes,'Forward',mainPath);
   dongCorrelFFTWOtime(char(registro{r}(2:end-4)),FOR_reshapeframes,'Forward',mainPath,0);
   dongCorrelFFTWOtime(char(registro{r}(2:end-4)),FOR_reshapeframes,'Forward',mainPath,1);
   
   disp('dongfourier-Backward');
   
   %dongCorrelFFT(char(registro{r}(2:end-4)),BACK_reshapeframes,'Backward',mainPath);
   dongCorrelFFTWOtime(char(registro{r}(2:end-4)),BACK_reshapeframes,'Backward',mainPath,0);
   dongCorrelFFTWOtime(char(registro{r}(2:end-4)),BACK_reshapeframes,'Backward',mainPath,1);
   
   clear ForwardMovieCRF;
   clear BackwardMovieCRF;
   clear FOR_reshapeframes;
   clear BACK_reshapeframes;
    
end


end

