function tmatrix = filtro(matrix,kernel)



lk = size(kernel,1);

mc = size(matrix,2);

line = zeros(1,mc);

for l=1:(lk-1)/2
    
   matrix = cat(1,matrix,line);
    
end

for l=1:(lk-1)/2
    
   matrix = cat(1,line,matrix);
    
end

mr = size(matrix,1);

column = zeros(1,mr);
column = column';

for l=1:(lk-1)/2
    
   matrix = cat(2,matrix,column);
    
end

for l=1:(lk-1)/2
    
   matrix = cat(2,column,matrix);
    
end

tmatrix = matrix;

ri = size(matrix,1);
ci = size(matrix,2);

for i=1:ri

    for j=1:ci
        
        imgij = 0;
        
        for h=-(lk-1)/2:(lk-1)/2

            for k=-(lk-1)/2:(lk-1)/2
                
                if ((i-h)>0 && (j-k)>0) && ((i-h)<(ri+1) && (j-k)<(ci+1))
                    
                    imgij = imgij + kernel((h+1+((lk-1)/2)),(k+1+((lk-1)/2)))*matrix(i-h,j-k);
                    
                end
                
            end
            
        end
        
        tmatrix(i,j) = imgij;
    
    end

end

end