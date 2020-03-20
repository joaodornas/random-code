function [Bx, By] = gradient(A)

x = size(A,1);
y = size(A,2);

for i=1:x
    
    for j=1:y
        
        if (i == x)
            
            k = 0;
            
        else
            
            k = 1;
            
        end
          
        Bx(i,j) = ( A(i + k,j) - A(i,j) ) / 2 ;
        
    end
    
end

for j=1:y
    
    for i=1:x
        
        if (j == y) 
            
            k = 0;
            
        else
            
            k = 1;
            
        end
          
        By(i,j) = ( A(i,j + k) - A(i,j) ) / 2 ;
        
    end
    
end


end

