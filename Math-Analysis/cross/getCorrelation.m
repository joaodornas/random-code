function correlation_norm = getCorrelation(A,B)



    sd_A = std(A);
    
    sd_B = std(B);
    
%     M_A = mean(A);
%     
%     M_B = mean(B);
   
%     A = A - M_A;
%     
%     B = B - M_B;

    dim_A = max(size(A,1),size(A,2));
    
    dim_B = max(size(B,1),size(B,2));
    
    dim = max(dim_A,dim_B);
    
    if dim_A > dim_B, B = [ B, zeros(1, (dim_A - dim_B) ) ]; end
    
    if dim_B > dim_A, A = [ A, zeros(1, (dim_B - dim_A) ) ]; end
    
    correlation = double(zeros(1, 2*dim - 2 ));
    
    correl = double(zeros(1,dim));
  
    for lag=(-dim+1):(dim-1)

               if lag ~= 0

                  for y=1:dim

                        if lag < 0

                            if (lag + y) < 0
                                
                                correl(y) = 0;
                                
                            elseif (lag + y) > 0

                                lag_v = lag + y;
    
                                correl(y) = A(y).*B(lag_v);

                            elseif (lag + y) == 0
                                
                                correl(y) = 0;
    
                            end

                        elseif lag > 0

                            if (lag + y) > dim
                                
                                correl(y) = 0;

                            elseif (lag + y) <= dim

                                lag_v = lag + y;

                                correl(y) = A(y).*B(lag_v);

                            end

                        end

                   end

                elseif lag == 0

                    correl = A.*B;
            
                end

            correlation(lag + dim) = mean(correl); 
        
            correl = [];

        end
    
    correlation_norm = correlation ./ (sd_A * sd_B) ;   

%     correlation_norm = abs(correlation_norm);

end

