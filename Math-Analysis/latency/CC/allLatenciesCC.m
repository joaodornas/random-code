function allLatenciesCC

registro = importdata('CC-Protocols.txt');

for r=56:73
    
   disp(strcat('Registro:',registro{r}(2:end-4))); 
    
   calcLatencyCC(char(registro{r}));
    
end

end

