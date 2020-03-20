function start = rateBefore(rate,h)
 
     FOUND = 'NO';
     
     for u=h:1

       if rate(u) == 0

            start = u;
            
            FOUND = 'YES';
            
            break;
            
       end

     end
     
      if (strcmp(FOUND,'YES') && ((start < (h - 3))))
         
         start = h - 3;
         
     end
     
     if strcmp(FOUND,'NO')
         
         start = h - 3;

         if (start < 0) || (start == 0)
            
            start = 1;

         end
         
     end

end