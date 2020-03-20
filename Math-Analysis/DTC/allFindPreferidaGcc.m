function allFindPreferidaGcc

registro = importdata('cc-44px-2sizes-g-Protocols.txt');

nConditions = 16;

angles = 22.5;

start_time = 0;

end_time = 3000;

for r=1:length(registro)
   
      results(r).name = char(registro{r}(2:end-4));
       
      [m, A, preferida, anti_preferida, rate, DI, response, sd] = getDirectionalSelectivity(registro{r},start_time,end_time);     
      
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
          print(f,'-depsc',strcat('/Users/joaodornas/Documents/_Research/_DATA/Center-Surround/DTC/',char(results(r).name),'-','vonMise-fit.eps'));

      else
          
          print(f,'-depsc',strcat('/Users/joaodornas/Documents/_Research/_DATA/Center-Surround/DTC/',char(results(r).name),'.eps'));

      end
      
      
       
end   


save(strcat('/Users/joaodornas/Documents/_Research/_DATA/Center-Surround/DTC/preferida'),'results');


end

