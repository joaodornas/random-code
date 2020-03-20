function end_ = rateAfter(rate,r,size_rate)

    FOUND = 'NO';
    
    for v=r:size_rate
    
       if rate(v) == 0

           end_ = v;
            
           FOUND = 'YES';
           
            break;

       end

    end
    
    if (strcmp(FOUND,'YES') && (end_ > (r + 3)))
          
         end_ = r + 3;
         
    end
     
     if strcmp(FOUND,'NO')
         
         end_ = r + min(4,length(r:size_rate)) - 1;
         
     end
               
end