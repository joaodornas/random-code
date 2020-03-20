function allFindPreferida

registro = importdata('memoryBackwardProtocols.txt');

nConditions = 16;

angles = 22.5;

for r=1:length(registro)
    
   [start_time, end_time, which_one] = DTC_forback(char(registro{r}));
   
   results(r).name = char(registro{r});
   
   if strcmp(which_one,'none')
       
       results(r).angle = NaN;
       
   else
       
      [m, A, preferida, anti_preferida, rate, DI, response, sd] = getDirectionalSelectivity(which_one,start_time,end_time);     
      
      results(r).preferida = preferida;
      
      results(r).anti_preferida = anti_preferida;
      
      results(r).A = A;
      
      results(r).m = m;
      
      results(r).rate = rate;
      
      results(r).DI = DI;
      
      results(r).response = response;
      
      results(r).sd = sd;
      
      if DI > 0.5
          
          for i=1:nConditions

                M(i) = vonMise(m,A,(i-1)*angles,preferida,anti_preferida,DI);

          end
          
      end
    
      f = figure;
      errorbar(1:nConditions,rate,2*sd);
      hold on;
      
      if DI > 0.5
          
          plot(1:16,M,'r');
          print(f,'-depsc',strcat('/Volumes/Data/DATA/Forward-Backward/DTC/',char(registro{r}(2:end-4)),'-',which_one(2:end-4),'-','vonMise-fit.eps'));

      else
          
          print(f,'-depsc',strcat('/Volumes/Data/DATA/Forward-Backward/DTC/',char(registro{r}(2:end-4)),'-',which_one(2:end-4),'.eps'));

      end
      
      
       
   end
    
end


save(strcat('/Volumes/Data/DATA/Forward-Backward/DTC/preferida'),'results');


end

