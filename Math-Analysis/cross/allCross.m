function allCross

registro = importdata('memoryBackwardProtocols.txt');

for r=1:length(registro)
    
    crossExp(char(registro{r}),3);
    
    %crossExp(char(registro{r}),2);
    
    close all;
    
end

end

