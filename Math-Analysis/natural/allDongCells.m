function allDongCells(start_registro,end_registro,OS)

registro = importdata('memoryBackwardProtocols.txt');

nFrames = 300;

save = 1;

if strcmp(OS,'Mac')
    
    mainPath = strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Dong/');
    
elseif strcmp(OS,'Win')
    
    mainPath = strcat('Z:\DATA\Forward-Backward\Power-Law\Dong\');
    
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
   
   for i=1:nFrames
    
        FOR_frame(:,:) = ForwardMovieCRF(:,:,i);
        
        BACK_frame(:,:) = BackwardMovieCRF(:,:,i);
   
        reshapeframe_FOR = reshape(FOR_frame.',[],1).';
        
        reshapeframe_BACK = reshape(BACK_frame.',[],1).';
    
        FOR_reshapeframes(:,i) = reshapeframe_FOR(:); 
        
        BACK_reshapeframes(:,i) = reshapeframe_BACK(:); 
    
   end
   
   dongCorrelFFT(char(registro{r}(2:end-4)),FOR_reshapeframes,'Forward',mainPath);
   
   dongCorrelFFT(char(registro{r}(2:end-4)),BACK_reshapeframes,'Backward',mainPath);
   
   clear ForwardMovieCRF;
   clear BackwardMovieCRF;
   clear FOR_reshapeframes;
   clear BACK_reshapeframes;
    
end


end

